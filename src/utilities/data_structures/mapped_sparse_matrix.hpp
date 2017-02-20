////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                           
//!                                                                             
//!                      \date Last Modified: 28/05/2014                        
//!                                                                             
//!	 \file
//!     This file contains a class to implement sparse array storage
//!     using a hash table scheme.
//!
//!     NOTES: 
//!
//!     The SPARSE MAPPED type storage scheme is a wrapper using the following
//!     classes. The class used can be set by a precompiler definition:       
//!
//!     1. std::unordered_map implements hash table data structure (c++ 11 
//!     is required). It is relatively slow to populate compared to
//!     other implementations and there is a high memory overhead.
//!     To add insult to injury, the memory clear function does not 
//!     reliably release memory resources back to the operating system. Avoid!
//!
//!     2. The sparsehash/sparse_hash_map library is highly memory efficient
//!     with a 2-bit per element overhead.   
//!
//!     3. The sparsehash/dense_hash_map library is highly speed efficient
//!     but with a much higher memory overhead compared to the sparse 
//!     implementation (it's still better than std::map however).
//!
//!     Both libraries 2 and 3 reliably release memory resources when the resize(0) 
//!     function is called. One disadvantage is that you need to reserve a map key
//!     to store deleted entries (also an empty map key in the dense map case)
//!     For convenience I choose the last possible map key 
//!     (e.g. the highest int or long int value) to be the deleted 
//!     entry key and (last possible value -1) for the empty map key. 
//!     Another issue is that the map is unordered, and it must be sorted if 
//!     converted to another format e.g. CRS.
//!
//!     http://research.neustar.biz/2011/11/27/big-memory-part-3-5-google-sparsehash/
//!                                                        
//!                    Copyright (C) 2014 Simon C Davenport
//!                                                                             
//!     This program is free software: you can redistribute it and/or modify
//!     it under the terms of the GNU General Public License as published by
//!     the Free Software Foundation, either version 3 of the License,
//!     or (at your option) any later version.
//!                                                                             
//!     This program is distributed in the hope that it will be useful, but
//!     WITHOUT ANY WARRANTY; without even the implied warranty of
//!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//!     General Public License for more details.
//!                                                                             
//!     You should have received a copy of the GNU General Public License
//!     along with this program. If not, see <http://www.gnu.org/licenses/>.
//!                                                                             
//////////////////////////////////////////////////////////////////////////////// 

#ifndef _MAPPED_SPARSE_MATRIX_HPP_INCLUDED_
#define _MAPPED_SPARSE_MATRIX_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////		 

#include "sparse_matrix.hpp"                //  For storageMethod_t
#include "../general/dcmplx_type_def.hpp"   //  For dcmplx type 
#include "../general/template_tools.hpp"    //  For is_same template function
#include "../general/cout_tools.hpp"        //  For manipulation of std::out
#include "../wrappers/mpi_wrapper.hpp"      //  Wrapper for mpi functionality
#include "../wrappers/murmur_hash_wrapper.hpp"//  For MurmurHasher128Wrapper
#include "../general/serialize.hpp"         //  For binary serialize functions

#include <limits>                           //  For finding the max 
                                            //  value of integer types

#if _SPEED_OPTIMIZED_MAP_
#include <sparsehash/dense_hash_map>    //  for google::dense_hash_map
#elif _MEMORY_OPTIMIZED_MAP_
#include <sparsehash/sparse_hash_map>   //  for google::sparse_hash_map
#else
#include <unordered_map>                //  STL unordered map if other
                                        //  implementations not available
                                        //  (requires c++11)
#endif

#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief Comparison operator to test if two CrsSparseMatrix data structures 
    //! are identical up to a given tolerence. Only classes with the same template
    //! arguments can be compared.
    //!
    //! Note that the comparison algorithm requires an in-place sort of the 
    //! underlying data and column arrays in the CrsSparseMatrix class. This
    //! does not change the data represented. 
    //!
    //////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool SameTest(
        MappedSparseMatrix<T>& lhs,     //!<    lhs matrix to compare
        MappedSparseMatrix<T>& rhs,     //!<    rhs matrix to compare
        const double tol)               //!<    comparison tolerence
    {
        if(
        (lhs.GetNbrNonZeros() != rhs.GetNbrNonZeros()) ||
        (lhs.GetNbrRows()     != rhs.GetNbrRows())     ||
        (lhs.GetNbrColumns()  != rhs.GetNbrColumns())
        )
        {
            return false;
        }
        
        if(utilities::is_same<dcmplx,T>::value)
        {
            for(auto it = lhs.m_map.begin(); it != lhs.m_map.end(); ++it)
            {
                uint32_t row;
                uint32_t col;
                
                utilities::Unpack2x32(it->first,row,col);

                if(abs(rhs.Value(row,col)-lhs.Value(row,col))>tol)  return false;
            }
        }
        else if(utilities::is_same<double,T>::value)
        {
            for(auto it = lhs.m_map.begin(); it != lhs.m_map.end(); ++it)
            {
                uint32_t row;
                uint32_t col;
                
                utilities::Unpack2x32(it->first,row,col);
                
                if(abs(rhs.Value(row,col)-lhs.Value(row,col))>tol)  return false;
            }
        }

        return true;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A class template to implement the mapped sparse matrix format
    //!
    //! This class uses the std::map container to represent a sparse matrix:
    //! The map associated a combined row/column index index i*nbrCols+j with
    //! a matrix element.
    //!
    //! Template parameter T is the matrix element type e.g. double of dcmplx
    //!
    //! NOTE: the default hash function for google::denst_hash_map or sparse_hash_map
    //! is std::tr1::hash<unsigned long>. The Murmur Hasher function is far superior
    //! 
    //////////////////////////////////////////////////////////////////////////////////
    
    template <typename T>
    class MappedSparseMatrix
    {
        private:
        
        //////////////////////////////////////////////////////////////////////////////////
        #if _SPEED_OPTIMIZED_MAP_
        
        google::dense_hash_map<uint64_t,T,MurmurHasher64Wrapper<uint64_t>> m_map;     
                                                //!<    An extremely fast hash map
                                                //!     implementation
                                               
        #elif _MEMORY_OPTIMIZED_MAP_
        
        google::sparse_hash_map<uint64_t,T,MurmurHasher64Wrapper<uint64_t>> m_map;    
                                                //!<    A memory efficient hash table
                                                //!     container to represent the matrix
        #else
        
        std::unordered_map<uint64_t,T> m_map;   //!<     Default to using the STL
                                                //!      unordered map
        #endif
        //////////////////////////////////////////////////////////////////////////////////

        int m_nbrRows;            //!<    Store number of matrix rows
        int m_nbrCols;            //!<    Store number of matrix columns

        static const uint64_t memoryQueryLimit = 0x100000000;
        //!<    Set 4GB max memory to store any given matrix before query
        
        template<typename R>
        friend class CrsSparseMatrix;
        
        template<typename R>
        friend bool SameTest(MappedSparseMatrix<R>& lhs,MappedSparseMatrix<R>& rhs,const double tol);

        //////////////////////////////////////////////////////////////////////////////////
        #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to return the default deleted map key
        //!
        //! Currently set to be the highest possible value 
        //! of the template integer type
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        uint64_t DeletedMapKey() const
        {
            return std::numeric_limits<uint64_t>::max();
        }
        
        #endif
        //////////////////////////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////////////////////////
        #if _SPEED_OPTIMIZED_MAP_ 
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to return the default empty map key
        //!
        //! Currently set to be the highest possible value 
        //! of the template integer type minus 1
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        uint64_t EmptyMapKey() const
        {
            return std::numeric_limits<uint64_t>::max()-1;
        }
        
        #endif
        //////////////////////////////////////////////////////////////////////////////////
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        public:
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Function to return estimated size of a map element
        //! (here we use that the 2-bit overhead may need to occupy a
        //! full word, hence the overhead is given a 4 byte word)
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        int SizeOfElement() const
        {
            return sizeof(uint64_t)+sizeof(T)+4;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Initializer from specified matrix dimensions
        //!
        //! map does not require memory pre-allocation
        //!
        //////////////////////////////////////////////////////////////////////////////////
    
        void Initialize(
            const int nbrRows,     //!<    Specify number of matrix rows
            const int nbrCols,     //!<    Specify number of matrix columns
            const long int nnz)    //!<    Specify estimated number of non-zeros
        {
            m_nbrRows = nbrRows;
            m_nbrCols = nbrCols;

            #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
            m_map.resize(nnz);
            #else
            m_map.reserve(nnz);
            #endif

            //  Print a warning message
        
            if(nnz*SizeOfElement()>memoryQueryLimit)
            {
                std::cerr<<"\n\tWARNING: POTENTIALLY EXCEEDING "<<memoryQueryLimit/(1024*1024*1024)<<" GB";
                std::cerr<<" MEMORY QUERY LIMIT! PRESS ANY KEY TO CONTINUE"<<std::endl;
                getchar();
            }
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//        
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Initialize from a dense matrix array
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void Initialize(
            T* denseFormat,             //!<    Dense version of matrix
            const int nbrRows,          //!<    Specify number of matrix rows
            const int nbrCols,          //!<    Specify number of matrix columns
            const double sparseZeroTol) //!<    Tolerance to assume matrix elements are zero
        {
            m_nbrRows = nbrRows;
            m_nbrCols = nbrCols;
        
            T* p_element = denseFormat;
    
            uint64_t dim = (uint64_t)nbrRows*(uint64_t)nbrCols;
    
            for(uint64_t i=0;i>dim;++i,++p_element)
            {
                if(abs((*p_element))>sparseZeroTol)
                {
                    m_map[i] = denseFormat[i];
                }
            }
            
            return;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Initialize from a CRS matrix. The data values must be the same type T
        //!
        //! NOTE 1: This initialization function will dealocate the original
        //! CrsSparseMatrix data structure in the process.
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void Initialize(
            CrsSparseMatrix<T>* other,  //!<    CRS format matrix to construct from
            const int nbrRows,          //!<    Specify number of matrix rows
            const int nbrCols)          //!<    Specify number of matrix columns
        {
            m_nbrRows = nbrRows;
            m_nbrCols = nbrCols;
        
            for(int i=0;i<other->m_rowStarts.size()-1;++i)
            {
                for(int j = other->m_rowStarts[i];j<other->m_rowStarts[i+1];++j)
                {
                    //m_map[i*m_nbrCols+other->m_cols[j]] =  other->m_data[j];
                    m_map[utilities::Pack2x32(i,other->m_cols[j])] =  other->m_data[j];
                }
            }

            for(int j = other->m_rowStarts.back();j<other->m_data.size();++j)
            {
                //m_map[(other->m_rowStarts.size()-1)*m_nbrCols+other->m_cols[j]] = other->m_data[j];
                m_map[utilities::Pack2x32(other->m_rowStarts.size()-1,other->m_cols[j])] = other->m_data[j];
            }

            other->Clear();
            
            return;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //!
        //! Default constructor
        //!
        
        MappedSparseMatrix()
        :
            m_nbrRows(0),
            m_nbrCols(0)
        {
            #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
            
            m_map.set_deleted_key(this->DeletedMapKey());
            
            #endif
            
            #if _SPEED_OPTIMIZED_MAP_
            
            m_map.set_empty_key(this->EmptyMapKey());
            
            #endif
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Copy constructor
        //!
        ////////////////////////////////////////////////////////////////////////////////// 
        
        MappedSparseMatrix(
            const MappedSparseMatrix* other)
            :
            m_map(other->m_map),
            m_nbrRows(other->m_nbrRows),
            m_nbrCols(other->m_nbrCols)
        {
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Destructor
        //!
        //! NOTE: with std::map (black-red binary tree) and std::unordered map (hash table)
        //! the clear function does not always free up memory resourced back to the OS!
        //!
        //! On the other hand the google::sparse_hash_map works fine as long as a deleted 
        //! key is defined (which is done in the constructors)
        //!
        ////////////////////////////////////////////////////////////////////////////////// 
        
        ~MappedSparseMatrix()
        {
            this->Clear();
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Clear internal memory allocations
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void Clear()
        {
            if(!m_map.empty())
            {
                m_map.clear();
                
                #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
                m_map.resize(0);
                #endif
            }
        }
       
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\// 
    
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Function to insert a matrix element into the map
        //! 
        //! \return The address of the matrix element, which can then be set with 
        //! an = operation e.g. Insert(i,j) = value
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        T& Insert(
            const int i,      //!<    row index
            const int j)      //!<    column index
        {
            //  Check the index is in range
            
            if(i<0 || i >= m_nbrRows || j<0 || j>= m_nbrCols)
            {
                std::cerr<<"\n\tERROR in MappedSparseMatrix with accessing element ["<<i<<","<<j<<"]"<<std::endl;
                std::cerr<<"\tINDEX OUT OF BOUNDS"<<std::endl;
                
                return m_map[0];
            }
            else
            {
                return m_map[utilities::Pack2x32(i,j)];
                //return m_map[i*m_nbrRows+j];
            }
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\// 
    
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Accessor function to read matrix element value (slow)
        //! 
        //////////////////////////////////////////////////////////////////////////////////
        
        const T Value(
            const int i,      //!<    row index
            const int j)      //!<    column index
            const
        {
            //  Check the index is in range
            
            if(i<0 || i >= m_nbrRows || j<0 || j>= m_nbrCols)
            {
                std::cerr<<"\n\tERROR in MappedSparseMatrix with accessing element ["<<i<<","<<j<<"]"<<std::endl;
                std::cerr<<"\tINDEX OUT OF BOUNDS"<<std::endl;
                
                return 0.0;
            }
            else
            {
                //auto it = m_map.find(i*m_nbrRows+j);
                auto it = m_map.find(utilities::Pack2x32(i,j));

                if(m_map.end() == it)
                {
                    //  Index not found in map, so return 0
                    
                    return 0.0;
                }
                else
                {
                    return it->second;
                }
            }
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\// 
    
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Define a function to convert from sparse-mapped to dense format
        //! 
        //////////////////////////////////////////////////////////////////////////////////

        void ConvertToArray(
            T* denseFormat)               //!<    Dense array output (with sufficient
                                          //!     memory allocated to store the full matrix)
            const
        {
            for(auto& x: m_map)
            {
                uint32_t row;
                uint32_t col;
                
                utilities::Unpack2x32(x.first,row,col);
                
                denseFormat[row*m_nbrRows+col] = x.second;
            }
        }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to remove any very small elements from the map
        //! to save on storage space
        //!
        ////////////////////////////////////////////////////////////////////////////////
        
        void Trim(
            const double sparseZeroTol)   //!<    Tolerance to assume matrix elements are zero
        {
            //  Iterate over the map and test for the zero condition

            auto it = m_map.begin();
            while(it != m_map.end())
            {
                if(utilities::is_same<T,dcmplx>::value)
                {
                    if(abs(it->second)<sparseZeroTol)
                    {
                        #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
                        
                        m_map.erase(++it);
                        
                        #else
                        
                        it = m_map.erase(it);
                        
                        #endif
                    }
                    else
                    {
                        ++it;
                    }
                }
                else if(utilities::is_same<T,double>::value)
                {
                    if(fabs(it->second)<sparseZeroTol)
                    {
                        #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
                    
                        m_map.erase(++it);
                    
                        #else
                        
                        it = m_map.erase(it);
                        
                        #endif
                    }
                    else
                    {
                        ++it;
                    }
                }
            }
        
            return;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Return the number of non-zero elements
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        long int GetNbrNonZeros() const
        {
            return m_map.size();
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the number of rows of the matrix
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        int GetNbrRows() const
        {
            return m_nbrRows;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the number of columns of the matrix
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        int GetNbrColumns() const
        {
            return m_nbrCols;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Print the non-zero matrix elements
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void Print() const
        {
            for(auto it = m_map.begin(); it != m_map.end(); ++it)
            {
                //utilities::cout.SecondaryOutput()<<"("<<floor(it->first/m_nbrCols)<<","<<it->first%m_nbrCols<<")\t";
                
                uint32_t row;
                uint32_t col;
                
                utilities::Unpack2x32(it->first,row,col);
                
                utilities::cout.SecondaryOutput()<<"\t\t("<<row<<","<<col<<")\t";

                utilities::cout.SecondaryOutput()<<it->second<<std::endl;
            }
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\// 

    };  //  End class MappedSparseMatrix
    
}   //  End namespace utilities

#endif
