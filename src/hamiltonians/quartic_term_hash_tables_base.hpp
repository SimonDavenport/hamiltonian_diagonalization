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

#ifndef _QUARTIC_TERM_HASH_TABLES_BASE_HPP_INCLUDED_
#define _QUARTIC_TERM_HASH_TABLES_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../utilities/general/orbital_and_state_defs.hpp"
#include "../utilities/data_structures/multi_key_hash.hpp"
#include "../utilities/wrappers/mpi_wrapper.hpp"
#include "../program_options/general_options.hpp"
#include <vector>

namespace diagonalization
{
    template<typename T>
    class QuarticTermHashTablesBase
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
        QuarticTermHashTablesBase()
        :
            m_kMax(0),
            m_fileHeader("# File contains term table data in hash table format")
        {}

        //!
        //! Destructor
        //!
        ~QuarticTermHashTablesBase(){}

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
        //! for a given k2,k3,k4. This allows us to pre-allocate a memory buffer
        //! to store the return values statically
        //!
        iSize_t GetMaxKCount() const
        {
            return m_kTable.GetMaxCount();
        }
        
        //!
        //! Get the dimension of the look-up table
        //!
        iSize_t GetDimension() const
        {
            return m_vTable.Size();
        }

        //!
        //! Get the k1 value for a given k2,k3,k4 with a 4 momentum conservation 
        //! law
        //!
        void GetK1(
            kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
            iSize_t& nbrK1,             //!<    Set the number of returned values             
            const kState_t k2,          //!<    Given k2 value
            const kState_t k3,          //!<    Given k3 value
            const kState_t k4)          //!<    Given k4 value
            const
        {
            m_kTable.Value(kRetrieveBuffer, nbrK1, utilities::Key(k2, k3, k4));
        }

        //!
        //! Return a specified coefficient of the quartic term
        //!
        T GetVkkkk(
            const kState_t k1,  //!<    Given k1 value
            const kState_t k2,  //!<    Given k2 value
            const kState_t k3,  //!<    Given k3 value
            const kState_t k4)  //!<    Given k4 value
            const
        {
            return m_vTable.Value(utilities::Key(k1, k2, k3, k4)); 
        }
        
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A template function to set the quartic hash table values from
        //! any other class container storing the quartic terms. 
        //!
        //! The function expects a class containing the methods:
        //!
        //! - GetK1(&k1,nbrK1,k2,k3,k4)
        //!
        //! Returns k1 for a given k2,k3,k4 momentum
        //!
        //! - GetVkkkk(k1,k2,k3,k4)
        //!
        //! Returns a single Vkkkk coefficient for a given k1,k2,k3,k4 set of momenta
        ////////////////////////////////////////////////////////////////////////////////
        template <class C>
        void SetFromArray(
            C* quarticArray,        //!<    Class container for quartic term array
            const kState_t filter)  //!<    Upper limit filter on k1
        {
            for(kState_t k2=0; k2<m_kMax; ++k2)
            {
                for(kState_t k3=0; k3<m_kMax; ++k3)
                {
                    for(kState_t k4=0; k4<m_kMax; ++k4)
                    {
                        kState_t k1;
                        iSize_t nbrK1;
                        quarticArray->GetK1(&k1, nbrK1, k2, k3, k4);
                        if(k1 <= filter)
                        {
                            m_kTable.Insert(k1, utilities::Key(k2, k3, k4));
                            m_vTable.Insert(utilities::Key(k1, k2, k3, k4)) = quarticArray->GetVkkkk(k1, k2, k3, k4);
                        }
                    }
                }
            }
        }
        
        //!
        //! Write quartic terms to a file
        //!
        void ToFile(
            const std::string fileName,     //!<    Name of file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class 
        {
            std::ofstream f_out;
            int nbrLabels = 4;
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
        //! Read quartic terms from a file
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
