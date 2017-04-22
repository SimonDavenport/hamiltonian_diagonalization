////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                  
//!                                                                             
//!	 \file
//!     This file contains the base class implementation for term tables
//!     using an underlying multi key hash table implementation
//!                                                  
//!                    Copyright (C) Simon C Davenport
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

#ifndef _QUADRATIC_TERM_HASH_TABLES_BASE_HPP_INCLUDED_
#define _QUADRATIC_TERM_HASH_TABLES_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../utilities/general/orbital_and_state_defs.hpp"
#include "../utilities/data_structures/multi_key_hash.hpp"
#include "../utilities/wrappers/mpi_wrapper.hpp"
#include "../program_options/general_options.hpp"
#include <vector>

namespace diagonalization
{
    template<typename T>
    class QuadraticTermHashTablesBase
    {
        protected:
        utilities::MultiHashMap<T> m_vTable;    //!< Look-up table data
        utilities::MultiHashMultiMap<kState_t> m_kTable;
                                            //!< Store a look-up table of k1 values 
                                            //!  satisfying a conservation law   
        kState_t m_kMax;                    //!<    Maximum index-value
        std::string m_fileHeader;           //!<    Text description at file head
        
        public:
        //!
        //! Default constructor
        //!
        QuadraticTermHashTablesBase()
        :
            m_kMax(0),
            m_fileHeader("# File contains term table data in hash table format")
        {}

        //!
        //! Destructor
        //!
        ~QuadraticTermHashTablesBase(){}

        //!
        //! Set the internal kMax variable
        //!
        void Initialize(
            const kState_t kMax)   //!<    No. k states
        {
            m_kMax = kMax;
        }

        //!
        //! Remove all the currently stored hash tables
        //!
        void Clear()
        {
            m_vTable.Clear();
            m_kTable.Clear();
        }
        
        //!
        //! Get the address of the quantum number hash table
        //!
        utilities::MultiHashMultiMap<kState_t>* GetKTable()
        {
            return &m_kTable;
        }   

        //!
        //! Get the address of the underlying hash table
        //!
        utilities::MultiHashMap<dcmplx>* GetVTable()
        {
            return &m_vTable;
        }

        //!
        //! Synchronize table contents with a given node
        //!
        void MpiSynchronize(
            const iSize_t syncId,               //!<    Node to synchronize with
            const utilities::MpiWrapper& mpi)   //!<    Instance of the MPI wrapper class
        {
            m_vTable.MpiSynchronize(syncId, mpi);
            m_kTable.MpiSynchronize(syncId, mpi);
        }

        //!
        //! Gather table contents onto given node
        //!
        void MpiGather(
            const iSize_t gatherId,             //!<    Node to gather map data
            const utilities::MpiWrapper& mpi)   //!<    Instance of the MPI wrapper class
        {
            m_vTable.MpiGather(gatherId, mpi);
            m_kTable.MpiGather(gatherId, mpi);
        }
        
        //!
        //! Get the highest number of k1 values that will be returned
        //! for a given k2. This allows us to pre-allocate a memory buffer
        //! to store the return values statically
        //!
        iSize_t GetMaxKCount() const
        {
            iSize_t maxCount = 1;
            for(kState_t k2=0; k2<m_kMax; ++k2)
            {
                maxCount = std::max(maxCount, (iSize_t)m_kTable.Count(k2));
            }
            return maxCount;
        }

        //!
        //! Get the k1 value for a given k2, with a 2 value conservation law
        //!
        void GetK1(
            kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
            iSize_t& nbrK1,             //!<    Set the number of returned values             
            const kState_t k2)          //!<    Given k2 value
            const
        {
            //  Get the corresponding list of k1 values
            m_kTable.Value(kRetrieveBuffer, nbrK1, k2);
        }

        //!
        //! Return a specified coefficient of the quadratic term
        //!
        T GetEkk(
            const kState_t k1,  //!<    k1 value
            const kState_t k2)  //!<    k2 value
            const
        {
            return m_vTable.Value(utilities::Key(k1, k2));
        }
        
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Incorporate an array containing quadratic coefficients
        //!
        //! The function expects a class containing the methods:
        //!
        //! - GetK1(&k1,nbrK1,k2)
        //!
        //! Returns k1 for a given k2 momentum
        //!
        //! - GetEkk(k1,k2)
        //!
        //! Returns a single Vkkkk coefficient for a given k1,k2 set of momenta
        ////////////////////////////////////////////////////////////////////////////////
        template <class C>
        void SetQuadraticFromArray(
            C* quadraticArray)    //!<    Class container for quadratic term array
        {
            for(kState_t k2=0; k2<m_kMax; ++k2)
            {
                kState_t k1;
                iSize_t nbrK1;
                quadraticArray->GetK1(&k1, nbrK1, k2);
                m_kTable.Insert(k1, k2);
                m_vTable.Insert(utilities::Key(k1, k2)) = quadraticArray->GetEkk(k1, k2);
            }
        }
        
        //!
        //! Write quadratic terms to a file
        //!
        void ToFile(
            const std::string fileName,     //!<    Name of file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class 
        {
            std::ofstream f_out;
            int nbrLabels = 2;
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                GenFileStream(f_out, fileName, format, mpi);
                if(!mpi.m_exitFlag)
                {
                    if(io::_BINARY_ == format)
                    {
                        iSize_t strSize = m_fileHeader.size();
                        f_out.write((char*)&(strSize), sizeof(iSize_t));
                        f_out.write(m_fileHeader.c_str(), strSize);
                    }
                    else if(io::_TEXT_ == format)
                    {
                        f_out << m_fileHeader << "\n";
                    }
                }
            }
            mpi.ExitFlagTest();
            m_kTable.ToFile(f_out, format, nbrLabels-1, mpi);
            m_vTable.ToFile(f_out, format, nbrLabels, mpi);
            if(f_out.is_open())
            {
                f_out.close();
            }
        }
        
        //!
        //! Read quadratic terms from a file
        //!
        void FromFile(
            const std::string fileName,     //!<    Name of file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class 
        {
            std::ifstream f_in;
            std::string fileHeader;
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                GenFileStream(f_in, fileName, format, mpi);
                if(!mpi.m_exitFlag)
                {
                    if(io::_BINARY_ == format)
                    {
                        iSize_t strSize;
                        f_in.read(reinterpret_cast<char*>(&strSize), sizeof(iSize_t));
                        fileHeader.resize(strSize);
                        f_in.read(&fileHeader[0], strSize);
                    }
                    else if(io::_TEXT_ == format)
                    {
                         getline(f_in, fileHeader);
                    }
                    if(m_fileHeader != fileHeader)
                    {
                        std::cerr << "\n\tFile read in error - wrong file header " << fileHeader << std::endl;
                        mpi.m_exitFlag = true;
                        f_in.close();
                    }
                }
            }
            mpi.ExitFlagTest();
            m_kTable.FromFile(f_in, format, mpi);
            m_vTable.FromFile(f_in, format, mpi);
            if(f_in.is_open())
            {
                f_in.close();
            }
            this->MpiSynchronize(0, mpi);
        }
    };
}   //  End namespace diagonalization 
#endif
