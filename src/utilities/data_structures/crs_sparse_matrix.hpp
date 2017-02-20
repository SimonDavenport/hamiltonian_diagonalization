////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 28/05/2014
//!                                                                             
//!	 \file
//!     This file contains a class to implement sparse array storage
//!     using a compressed row storage scheme.
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

#ifndef _CRS_SPARSE_MATRIX_HPP_INCLUDED_
#define _CRS_SPARSE_MATRIX_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "sparse_matrix.hpp"                //  For storageMethod_t
#include "../general/dcmplx_type_def.hpp"   //  For dcmplx type 
#include "../general/template_tools.hpp"    //  For is_same template function
#include "../general/cout_tools.hpp"        //  For manipulation of std::out
#include "../wrappers/mpi_wrapper.hpp"      //  Wrapper for mpi functionality
#include "../algorithms/quick_sort.hpp"     //  For QuickSort function
#include "../general/serialize.hpp"         //  For binary serialize functions

#include <vector>                           //  For std::vector
#include <algorithm>                        //  For std::remove_if, 
                                            //  for binary_search and sort
#include <limits>                           //  For finding the max 
                                            //  value of integer types
#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{
    //////////////////////////////////////////////////////////////////////////////////
    //!	\brief Used for sorting paired lists
    //!	
    //////////////////////////////////////////////////////////////////////////////////

    template<typename U, typename V>
    struct pair
    {
        U first;
        V second;
        
        bool operator<(const pair<U,V>& val) const {return first < val.first;}
    };
    
    
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
        CrsSparseMatrix<T>& lhs,    //!<    lhs matrix to compare
        CrsSparseMatrix<T>& rhs,    //!<    rhs matrix to compare
        const double tol)           //!<    comparison tolerence
    {
        if(
        (lhs.GetNbrNonZeros() != rhs.GetNbrNonZeros()) ||
        (lhs.m_rowStarts      != rhs.m_rowStarts)
        )
        {
            return false;
        }
    
        //  In CRS format the row data is always stored in order

        for(int i=0;i<lhs.GetNbrRows();++i)
        {
            //  For each row, we need to ensure that the column indices are 
            //  in the same order (since this is not necessarily guaranteed)
            //  So run a quick sort to rearrange the columns and their 
            //  corresponding data entries
            
            int startCol = lhs.m_rowStarts[i];
            int endCol   = lhs.m_rowStarts[i+1];

            utilities::QuickSort<int,T,_ASCENDING_ORDER_>(&lhs.m_cols[startCol],&lhs.m_data[startCol],endCol-startCol);
            
            utilities::QuickSort<int,T,_ASCENDING_ORDER_>(&rhs.m_cols[startCol],&rhs.m_data[startCol],endCol-startCol);
        }

        //  Now test for equality of the underlying data structures
        
        if(lhs.m_cols != rhs.m_cols)
        {
            return false;
        }

        if(utilities::is_same<dcmplx,T>::value)
        {
            for(int i=0;i<lhs.GetNbrNonZeros();++i)
            {
                if(abs(lhs.m_data[i]-rhs.m_data[i])>tol)    return false;
            }
        }
        else if(utilities::is_same<double,T>::value)
        {
            for(int i=0;i<lhs.GetNbrNonZeros();++i)
            {
                if(fabs(lhs.m_data[i]-rhs.m_data[i])>tol)   return false;
            }
        }
        
        return true;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A template class to implement the compressed row storage 
    //! sparse matrix format
    //!
    //! This class uses std::vector for dynamic memory allocation
    //!
    //! The T argument is a matrix element type e.g. double or double complex
    //! 
    //////////////////////////////////////////////////////////////////////////////////
    
    template<typename T>
    class CrsSparseMatrix
    {
        private:
        
        std::vector<int> m_cols;        //!<    List of column indices, of dimension nbrNonZeros
        std::vector<int> m_rowStarts;   //!<    List containing the position of the start of each 
                                        //!     row in the list m_cols, of dimension nbrRows+1
                                        //!     The last value contains the number
                                        //!     of non-zero elements for simplified iteration
        std::vector<T> m_data;          //!<    List of data associated with each entry, of 
                                        //!     dimension nbrNonZeros

        static const uint64_t memoryQueryLimit = 0x100000000;
        //!<    Set 4GB max memory to store any given matrix before query
        
        bool m_preparedForParallelMultiply; //!<    A flag set to true only once the
                                            //!     PrepareForParallelMultiplication
                                            //!     function has been called
    
        template<typename R>
        friend class MappedSparseMatrix;

        template<typename R>
        friend bool SameTest(CrsSparseMatrix<R>& lhs,CrsSparseMatrix<R>& rhs,const double tol);

        public:
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief  Default constructor
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        CrsSparseMatrix()
        :
            m_preparedForParallelMultiply(false)
        {}
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Initialize with preallocation specified
        //!
        //! Note: when used in a parallel implementation nbrRows separately
        //! specifies the number of rows stored on each node
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void Initialize(
            const int nbrRows,      //!<    Specify number of rows
            const int nbrCols,      //!<    Specify number of columns
            const long int nnz)     //!<    Specify estimated number of non-zeros
                                    //!     (for dynamic memory reservation)
        {
            //  Reserve space in memory to store the array
        
            if((nnz*(sizeof(int)+sizeof(T))+(nbrRows+1)*sizeof(int))>memoryQueryLimit)
            {
                std::cerr<<"\n\tWARNING: EXCEEDING "<<memoryQueryLimit/(1024*1024*1024)<<" GB";
                std::cerr<<" MEMORY QUERY LIMIT! PRESS ANY KEY TO CONTINUE"<<std::endl;
                getchar();
            }
        
            //  Reserve an estimated amount of memory space to store the 
            //  column and data values
        
            m_cols.reserve(nnz);
            m_data.reserve(nnz);
            
            //  Generate an array of row start "pointers" initialized to zero
            
            m_rowStarts.resize(nbrRows+1);    
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Initialize from MappedSparseMatrix. The data values must be the same 
        //! type T
        //!
        //! NOTE 1: This initialization function will dealocate the original
        //! MappedSparseMatrix data structure in the process. 
        //!
        //! NOTE 2: the data in the hash table implementation is assumed to be unordered
        //! in general (it would be more efficient if it were ordered). This function
        //! performs a sorting operation on the map elements before converting to CRS format 
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void Initialize(    
            MappedSparseMatrix<T>* other)         //!<    Mapped matrix data structure
        {
            //  First call the std::sort algorithm on the unordered list
            //  of non-zero matrix elements
            
            //  We need to extract the map data as an array of pairs
            
            //  Allocate memory to store the ENTIRE map
            
            int nnz = other->m_map.size();

            pair<uint64_t,T>* sortingBuffer = new pair<uint64_t,T>[nnz];
            
            pair<uint64_t,T>* p_buffer = sortingBuffer;
            
            //  Set the pair values 
            
            for(auto& x: other->m_map)
            {
                p_buffer->first = x.first;
                p_buffer->second = x.second;

                ++p_buffer;
            }

            //  Deallocate memory associated with the original map 
            
            other->m_map.clear();
            
            #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
                other->m_map.resize(0);
            #endif

            //  Apply a std::sort algorithm on the array of pairs
            //  (sorting in order of the map key values)
            
            std::sort(sortingBuffer,sortingBuffer+nnz);
            
            //  Now convert the sorted arrays into CRS format
            
            m_cols.reserve(nnz);
            m_data.reserve(nnz);

            //  Generate an array of row start "pointers" initialized to zero

            m_rowStarts.resize(other->m_nbrRows+1);
            m_rowStarts.back() = nnz;

            int rowCounter = 0;
            int currRow = 0;
            p_buffer = sortingBuffer;

            for(int i=0;i<nnz;++i,++p_buffer)
            {
                uint32_t row;
                uint32_t col;

                utilities::Unpack2x32(p_buffer->first,row,col);

                //m_cols.push_back(temp%other->m_nbrCols);
                m_cols.push_back(col);
                m_data.push_back(p_buffer->second);

                //if(floor((double)temp/other->m_nbrCols) >= currRow)
                if(row >= currRow)
                {
                    //  Increment to the actual row value
                
                    //int actualRow = floor((double)temp/other->m_nbrCols);
                    //int actualRow = row;
                
                    while(currRow<=row)
                    {
                        //  Add to list of row starts
                        m_rowStarts[currRow] = rowCounter;
                        ++currRow;
                    }
                }

                ++rowCounter;
            }
            
            //  Set the remaining row starts to the final counter value
            
            while(currRow<other->m_nbrRows)
            {
                m_rowStarts[currRow] = rowCounter;
                
                ++currRow;
            }
            
            //  Deallocate the sorting buffer
            
            delete[] sortingBuffer;

            return;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //!  \brief Initialize from dense matrix
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void Initialize(
            T* array,                       //!<    Dense array
            const int  nbrRows,             //!<    Specify number of rows
            const int  nbrCols,             //!<    Specify number of columns
            const double sparseZeroTol)     //!<    Tolerance to assume matrix 
                                            //!     elements are zero
        {
            
            //  Preallocate the exact memory required for the conversion
            
            T* p_element = array;
            long int counter = 0;
        
            for(int i=0;i<nbrRows;++i)
            {
                for(int j=0;j<nbrCols;++j,++p_element)
                {
                    if(abs((*p_element))>sparseZeroTol)
                    {
                        ++counter;
                    }
                }
            }

            m_cols.reserve(counter);
            m_data.reserve(counter);
            
            m_rowStarts.resize(nbrRows+1);
            
            m_rowStarts.back() = counter;
            
            //  Iterate over the dense array and convert to CRS format
        
            int currRow = 0;
            counter   = 0;
            p_element = array;
            
            for(int i=0;i<nbrRows;++i)
            {
                for(int j=0;j<nbrCols;++j,++p_element)
                {
                    if(abs((*p_element))>sparseZeroTol)
                    {
                        m_cols.push_back(j);      
                        m_data.push_back(*p_element);
                        
                        if(i >= currRow)
                        {
                            //  Add to list of row starts
                            
                            m_rowStarts[currRow] = counter;
                            
                            ++currRow;
                        }
                        
                        ++counter;
                    }
                }
            }
            
            this->Print();
        }    
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Copy constructor
        //!
        ////////////////////////////////////////////////////////////////////////////////// 
        
        CrsSparseMatrix(
            const CrsSparseMatrix* other)
            :
            m_cols(other->m_cols),
            m_rowStarts(other->m_rowStarts),
            m_data(other->m_data),
            m_preparedForParallelMultiply(other->m_preparedForParallelMultiply)
        {
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Destructor
        //!
        //! In order to ensure that memory space is properly released, one 
        //! technique is to perform a swap operation with an empty local vector
        //!
        ////////////////////////////////////////////////////////////////////////////////// 
        
        ~CrsSparseMatrix()
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
            m_cols.clear();
            m_rowStarts.clear();
            m_data.clear();
            m_preparedForParallelMultiply = false;
        }
       
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \breif Accessor function to read matrix element value
        //!
        //! Use as e.g.: element = object.value(i,j); or std::cout<<object.value(i,j); etc.
        //!
        //! \return The value of the matrix element i,j (or 0.0 if the
        //! element could not be found)
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        const T Value(
            const int i,      //!<    row index
            const int j)      //!<    column index
            const
        {
            //  Since the array is unordered in general, we cannot use
            //  a fast binary search to find the matrix element
            //  (use the mapped format for fast indexing)
        
            if(i<0 || i >= m_rowStarts.size()-1)
            {
                std::cerr<<"\n\tERROR in CrsSparseMatrix with accessing element ["<<i<<","<<j<<"]"<<std::endl;
                std::cerr<<"\tINDEX OUT OF BOUNDS"<<std::endl;
                
                return 0.0;
            }
            else
            {
                for(int t = m_rowStarts[i];t<m_rowStarts[i+1];++t)
                {
                    if(m_cols[t] == j)  return m_data[t];
                }
            }
        
            //  If it was not found, then the element is zero
        
            return 0.0;
        }
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////   
        //! \brief Accessor function to insert or update matrix element value
        //!
        //! Use as: object.insert(i,j) = element; OR object.insert(i,j) += element; etc.
        //!
        //! \return Address of the matrix element i,j OR
        //! if that could not be found, a new element is added at that
        //! position and the address of the new element is returned. 
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        T& Insert(
            const int i,      //!<    row index
            const int j)      //!<    column index
        {
            //  Since the array is unordered in general, we cannot use
            //  a fast binary search to find the matrix element
            //  (use the mapped format for fast indexing)

            if(i<0 || i >= m_rowStarts.size()-1)
            {
                std::cerr<<"\n\tERROR in CrsSparseMatrix with accessing element ["<<i<<","<<j<<"]"<<std::endl;
                std::cerr<<"\tROW INDEX OUT OF BOUNDS"<<std::endl;
                
                return m_data[0];
            }
            else
            {
                for(int t = m_rowStarts[i];t<m_rowStarts[i+1];++t)
                {
                    if(m_cols[t] == j)  return m_data[t];
                }
                
                //  If it was not found, then insert a new element into 
                //  the list and return its value

                m_cols.insert(m_cols.begin()+m_rowStarts[i+1],j);
                m_data.insert(m_data.begin()+m_rowStarts[i+1],0.0);

                //  Bump up the remaining row starts
                
                for(int r = i+1;r<m_rowStarts.size();++r)
                {
                    ++m_rowStarts[r];
                }

                return m_data[m_rowStarts[i+1]-1];
            }
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////   
        //! \brief MPI gather function
        //!
        //! Combine instances of this class together onto a single instance on
        //! the master node. Memory allocations on remaining nodes are then cleared.
        //!
        //! NOTE: the row indices on each node are assumed to be independent. 
        //! So e.g. rows 1-100 on node 0 + rows 1-100 on node 1 become rows
        //! 1-200 on node 0.
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void MpiGather(
            const int gatherNode,   //!<    Node to gather data onto
            const MpiWrapper& mpi)  //!<    Instance of the MPI wrapper class
        {
            //  First we determine how many rows there will be in the combined matrix
            //  With the CRS Format in parallel, it is assumed that some number 
            //  of matrix rows are stored on each node in order

            int totalNbrRows = 0;
            int nbrRows = m_rowStarts.size()-1;

            MPI_Reduce(&nbrRows,&totalNbrRows,1,mpi.GetType<int>(),MPI_SUM,gatherNode,mpi.m_comm);

            //  Second we need the total number of non-zero elements and the 
            //  cumulative number stored on each node in order

            int nbrNonZeros = m_data.size();
            int totalNbrNonZero = 0;
            int allNbrNonZero[mpi.m_nbrProcs];
            int cumulativeNbrNonZero = 0;

            MPI_Reduce(&nbrNonZeros,&totalNbrNonZero,1,mpi.GetType<int>(),MPI_SUM,gatherNode,mpi.m_comm);

            MPI_Allgather(&nbrNonZeros,1,mpi.GetType<int>(),&allNbrNonZero,1,mpi.GetType<int>(),mpi.m_comm);

            for(int k=0;k<mpi.m_id;++k)
            {
                cumulativeNbrNonZero += allNbrNonZero[k];
            }

            //  Shift the row starts by the cumulative number of non-zeros
            //  on preceding nodes
            
            for(int i=0;i<nbrRows;++i)
            {
                m_rowStarts[i] += cumulativeNbrNonZero;
            }

            MPI_Barrier(mpi.m_comm);

            //  Next we combine the index and data arrays from each node
            //  It is assumed that each node contains a set of matrix rows
            //  in the correct order, and no row sorting takes place in the
            //  gathered matrix
            
            //  Pre-allocate the required memory to store the whole matrix
            //  on the master node

            if(gatherNode == mpi.m_id)
            {
                utilities::cout.DebuggingInfo()<<"\n\t- ON GATHER NODE: REALLOCATING (ESTIMATED) "<<((totalNbrRows*sizeof(int)+totalNbrNonZero*(sizeof(T)+sizeof(int)))/(1024.0*1024.0));
                utilities::cout.DebuggingInfo()<<" MB TO STORE ENTIRE SPARSE MATRIX"<<std::endl;

                m_rowStarts.resize(totalNbrRows+1);
                m_cols.resize(totalNbrNonZero);
                m_data.resize(totalNbrNonZero);
                
                m_rowStarts.back() = totalNbrNonZero;
            }

            MPI_Status status;

            mpi.Gather<int>(m_rowStarts.data(),nbrRows,
            m_rowStarts.data(),totalNbrRows,gatherNode,mpi.m_comm,status);
            
            mpi.Gather<int>(m_cols.data(),nbrNonZeros,
            m_cols.data(),totalNbrNonZero,gatherNode,mpi.m_comm,status);

            mpi.Gather<T>(m_data.data(),nbrNonZeros,
            m_data.data(),totalNbrNonZero,gatherNode,mpi.m_comm,status);

            MPI_Barrier(mpi.m_comm);
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////   
        //! \brief A function to resize the row starts array such that the stored
        //! matrix now contains an additional number of empty rows after the 
        //! populated part - assuming that the matrix is stored 
        //! in parallel (i.e. with some number of rows stored on each node).
        //! This turns it into a square matrix, which is required for the matrix
        //! vector multiplication routine. 
        //!
        //! The picture to have in mind is as follows. We solve Y = A X, where
        //! the memory distribution of the vectors/matrix is as:
        //!
        //!         +------------+
        //! Y0      | A0         |     X0   }   dataDistribution[0]
        //!         |   +--------+
        //! Y1      |   | A1     |     X1   }   dataDistribution[1]
        //!         |   |   +----+
        //! Y2      |   |   | A2 |     X2   }   dataDistribution [2]
        //!
        //! Where the first dataDistribution[0] rows of the matrix are stored on 
        //! node 0, and so on. Note that the Hermitian structure of the matrix 
        //! means that the first n rows also give the first n columns.
        //! 
        //! The matrix A0 can be thought of as a NxN matrix (N the total dimension)
        //! with zeros where A1, A2 are labelled. A1 is a (N-dataDistribution[0])x 
        //! (N-dataDistribution[0]) matrix with zeros where A2 is labelled, and so on.
        //! Appropriate restructuring is done by PrepareForParallelMultiplication.
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void PrepareForParallelMultiplication(
            std::vector<int>* dataDistribution, //!<    A vector to store the number
                                                //!     of rows on each node
            std::vector<MPI_Comm>* commGroups,  //!<    List of processor groups
                                                //!     required for minimal
                                                //!     vector MPI communication
            std::vector<T>* inputVectorBuffer,  //!<    Buffer for MPI distribution of X vector
            std::vector<T>* outputVectorBuffer, //!<    Buffer for MPI accumulation of Y vector
            const MpiWrapper& mpi)              //!<    Instance of the MPI wrapper class
        {
            //  First we determine how many rows are stored on each node, the
            //  the cumulative number stored and the total number stored. This
            //  information is required on all nodes.

            if(!m_preparedForParallelMultiply)
            {
                int cumulativeRows = 0;
                int totalRows;
                int nbrRows = m_rowStarts.size()-1;
                
                dataDistribution->resize(mpi.m_nbrProcs);

                MPI_Allreduce(&nbrRows,&totalRows,1,mpi.GetType<int>(),MPI_SUM,mpi.m_comm);
                
                MPI_Allgather(&nbrRows,1,mpi.GetType<int>(),&(*dataDistribution)[0],1,mpi.GetType<int>(),mpi.m_comm);
     
                for(int k=0;k<mpi.m_id;++k)
                {
                    cumulativeRows += (*dataDistribution)[k];
                }

                int newSize = totalRows-cumulativeRows;

                //  insert an array for the empty rows below the populated part
                
                int nnz = m_rowStarts.back();
      
                m_rowStarts.resize(newSize+1,nnz);
            
                //  Update the m_cols array to shrink the matrix size stored on each 
                //  node

                for(int i=0;i<nnz;++i)
                {
                    m_cols[i] -= cumulativeRows;
                }
                
                //  Define a list of processor groups that will allow
                //  minimal MPI communication during vector multiplication
                
                for(int i = 0;i<mpi.m_nbrProcs;++i)
                { 
                    //  Broadcast the contents of the buffer to the preceding i nodes 
                    //  only, by splitting up the MPI communicator into 2 groups -
                    //  one for the preceding nodes, one for the rest
                    
                    MPI_Comm  minimalComm;

                    int whichGroup = (i>=mpi.m_id ? 0 : 1);

                    MPI_Comm_split(mpi.m_comm,whichGroup,mpi.m_id,&minimalComm);

                    commGroups->push_back(minimalComm);
                }
                
                inputVectorBuffer->resize(newSize);
                outputVectorBuffer->resize(newSize);
                
                m_preparedForParallelMultiply = true;
            }
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
      
        //////////////////////////////////////////////////////////////////////////////////   
        //! \brief Print the non-zero elements of the sparse matrix
        //!
        //////////////////////////////////////////////////////////////////////////////////   
        
        void Print() const
        {   
            for(int i=0;i<m_rowStarts.size()-1;++i)
            {
                utilities::cout.SecondaryOutput()<<"\n\t";
                
                for(int j = m_rowStarts[i];j<m_rowStarts[i+1];++j)
                {
                    utilities::cout.SecondaryOutput()<<"\t\t("<<i<<","<<m_cols[j]<<")\t";
                    
                    utilities::cout.SecondaryOutput()<<m_data[j]<<"\t";
                }
            }

            return;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
#if 0        
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief	A function to evaluate the Frobenius norm of a sparse matrix
        //!
        //! \return  \[ \sqrt( \sum_{i,j} |A_ij|^2 ) \]
        //!
        ////////////////////////////////////////////////////////////////////////////////
        
        double FrobeniusNorm() const
        {
            double norm = 0.0;
    
            for(int i=0;i<m_rowStarts.size()-1;++i)
            {
                for(int j = m_rowStarts[i];j<m_rowStarts[i+1];++j)
                {
                    norm += m_data[j]*std::conj(m_data[j]);
                }
            }
            
            return sqrt(norm);
        }
#endif       
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the number of non-zero elements
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        long int GetNbrNonZeros() const
        {
            return m_data.size();
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the number of rows of the matrix
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        int GetNbrRows() const
        {
            return m_rowStarts.size()-1;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the const address of the start of the m_rowStarts array
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        const int* GetRowStarts() const
        {
            return &m_rowStarts[0];
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the address of the start of the m_rowStarts array
        //!
        //////////////////////////////////////////////////////////////////////////////////
                
        int* GetRowStarts()
        {
            return &m_rowStarts[0];
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the const address of the start of the m_cols array
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        const int* GetCols() const
        {
            return &m_cols[0];
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the address of the start of the m_cols array
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        int* GetCols()
        {
            return &m_cols[0];
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the const address of the start of the m_data array
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        const T* GetData() const
        {
            return &m_data[0];
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
             
        //////////////////////////////////////////////////////////////////////////////////
        //!\brief Return the address of the start of the m_data array
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        T* GetData()
        {
            return &m_data[0];
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief  Return dense matrix
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void ConvertToArray(
            T* denseMatrix,         //!<    Address of dense matrix to be returned
                                    //!     (allocated to be the appropriate size)
            const int nbrColumns)   //!<    Number of columns in the array
            const
        {
            for(int i=0;i<m_rowStarts.size()-1;++i)
            {
                for(int j = m_rowStarts[i];j<m_rowStarts[i+1];++j)
                {
                    denseMatrix[i*nbrColumns+m_cols[j]] =  m_data[j];
                }
            }
            
            return;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief  Place matrix data in a file
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void DataToFile(std::ofstream& f_out)
            const
        {
            f_out<<m_data.size()<<"\n";
            f_out<<m_rowStarts.size()-1<<"\n";

            for(int i=0;i<m_rowStarts.size()-1;++i)
            {
                for(int j = m_rowStarts[i];j<m_rowStarts[i+1];++j)
                {
                    if(utilities::is_same<double,T>::value)
                    {
                        f_out<<i<<" "<<m_cols[j]<<" "<<m_data[j]<<"\n";
                    }
                    else if(utilities::is_same<dcmplx,T>::value)
                    {
                        f_out<<i<<" "<<m_cols[j]<<" "<<std::real(m_data[j])<<" "<<std::imag(m_data[j])<<"\n";
                    }
                }
            }
            
            return;
        }
      
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to remove any very small elements from the matrix
        //! to save on storage space (slow)
        //!
        //! This function uses c++11 methods - compile with the --c++11 option 
        //!
        ////////////////////////////////////////////////////////////////////////////////
     
        void Trim(
            const double sparseZeroTol)   //!<    Tolerance to assume matrix elements are zero
        {
            //  Iterate over matrix elements and test for zero values

            //  First update the column index list and list of row starts

            //  Declare variables to pass to a c++11 "lambda function"

            std::vector<int>  newRowStarts(m_rowStarts);
            auto data_it = m_data.begin();
            auto row_it  = m_rowStarts.begin();
            auto newRow_it  = newRowStarts.begin();
            int nbrRows  = m_rowStarts.size()-1;
            int counter  = 0;

            auto it = std::remove_if(m_cols.begin(),m_cols.end(),
            
            //  Declare lambda function [variables](input){function body}
            [sparseZeroTol,&data_it,&counter,&row_it,&newRow_it,nbrRows](int x) -> bool
            {
                bool tmp = false;
            
                if(utilities::is_same<T,dcmplx>::value)
                {
                    tmp = abs(*data_it) < sparseZeroTol;
                }
                else if(utilities::is_same<T,double>::value)
                {
                    tmp = fabs(*data_it) < sparseZeroTol;
                }
                
                if(tmp)
                {
                    for(int r=0;r<nbrRows-1;++r)
                    {
                        if(counter>=row_it[r] && counter<row_it[r+1])
                        {
                            for(int k=r+1;k<nbrRows;++k)
                            {
                                newRow_it[k]--;
                            }
                            
                            break;
                        }
                    }
                }

                ++data_it;
                ++counter;
                
                return tmp;
            }
            );  //  Complete the remove_if call
            
            //  Update the list of row start positions
            
            m_rowStarts = newRowStarts;
            
            //  Clear the temporary list
            
            newRowStarts.clear();
            
            m_cols.erase(it,m_cols.end());

            //  Also update the data values

            auto it1 = std::remove_if(m_data.begin(),m_data.end(),
            [sparseZeroTol](dcmplx x){return abs(x) < sparseZeroTol;});

            m_data.erase(it1,m_data.end());

            //  Finally, shrink the containers to exactly fit the remaining data
            
            m_data.shrink_to_fit();
            m_cols.shrink_to_fit();
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
  
    };  //  End class CrsSparseMatrix
    
}   // End namespace utilities

#endif
