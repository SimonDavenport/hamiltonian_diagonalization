////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store data structures describing 
//!     quadratic terms in a Hamiltonian
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

#ifndef _QUADRATIC_TERM_TABLES_BASE_HPP_INCLUDED_
#define _QUADRATIC_TERM_TABLES_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../utilities/general/orbital_and_state_defs.hpp"
#include "../utilities/wrappers/io_wrapper.hpp"
#include "../utilities/wrappers/mpi_wrapper.hpp"
#include "../program_options/general_options.hpp"
#include <vector>
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    template<typename T>
    class QuadraticTermTablesBase
    {
        protected:
        std::vector<T> m_vTable;            //!<    Look-up table data 
        std::vector<kState_t> m_kTable;     //!<    Store a look-up table of k1 values 
                                            //!     satisfying a conservation law
        kState_t m_kMax;                    //!<    Maximum k-value
        std::string m_fileHeader;           //!<    Text description at file head
        
        //!
        //! Compute the dimension of the look-up table
        //!
        iSize_t CalculateDim(
            const kState_t kMax)    //!<    Number of k states
            const
        {
            return kMax;
        }
        
        public:
        //!
        //! Default constructor
        //!
        QuadraticTermTablesBase()
        :
            m_kMax(0),
            m_fileHeader("# File contains term table data in array format")
        {}

        //!
        //! Destructor
        //!
        ~QuadraticTermTablesBase()
        {
            this->Clear();
        }

        //!
        //! Initialize memory allocations
        //!
        void Initialize(
            const kState_t kMax,                //!<    No. k states
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class 
        {
            m_kMax = kMax;
            iSize_t dim = this->CalculateDim(m_kMax);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            { 
                utilities::cout.AdditionalInfo()<<"\n\t- INITIALIZING QUADRATIC TERM TABLE "<<std::endl;
                utilities::cout.DebuggingInfo()<<"\n\t- ALLOCATING "<<2*sizeof(T)*dim/(1024.0*1024.0)
                         <<" mb to store look-up table of "<<dim<<" terms."<<std::endl;
            }
            m_kTable.resize(dim);
            m_vTable.resize(dim);
        }
        
        //!
        //! Reset the memory allocations.
        //!
        void Clear()
        {
            m_kTable.clear();
            m_vTable.clear();
        }
        
        //!
        //! Only allow for a single k value to be returned
        //!
        iSize_t GetMaxKCount() const
        {
            return 1;
        }

        //!
        //! Get the dimension of the look-up table
        //!
        iSize_t GetDimension() const
        {
            return m_vTable.size();
        }

        //!
        //! Return the address of the class momentum table so that 
        //! it can be set externally.
        //!
        std::vector<kState_t>* GetKTable()
        {
            return &m_kTable;
        }

        //!
        //! Return the address of the class quadratic term array so 
        //! that it can be set externally.
        //!
        std::vector<T>* GetVTable()
        {
            return &m_vTable;
        }

        //!
        //! Given an overall coefficient g multiplying all of the 
        //! quartic terms, this function incorporates the coefficient into 
        //! all table values.
        //!
        void SetCoefficient(
            const T coefficient)   //!<    Coefficient value
        {
            for(auto& it : m_vTable)
            {
                it *= coefficient;
            }
        }

        //!
        //! Synchronize table contents with a given node
        //!
        void MpiSynchronize(
            const iSize_t syncId,               //!<    Node to synchronize with
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            mpi.Sync(&m_vTable, syncId);
            mpi.Sync(&m_kTable, syncId);
        }
        
        //!
        //! Find the k1 value corresponding to a given k2.
        //!
        void GetK1(
            kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
            iSize_t& nbrK1,             //!<    Set the number of returned values 
            const kState_t k2)          //!<    k2 index
            const
        {
            nbrK1 = 1;
            //  In this case the quadratic terms are diagonal
            kRetrieveBuffer[0] =  k2;
        }
        
        //!
        //! Return a quadratic term coefficient for a given k1,k2
        //!
        T GetEkk(
            const kState_t k1,         //!<    k1 index
            const kState_t k2)         //!<    k2 index
            const
        {
            return this->m_vTable[k1];
        }
        
        //!
        //! Write quadratic terms to a file
        //!
        void ToFile(
            const std::string fileName,     //!<    Name of file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class 
            const
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                int nbrLabels = 2;
                std::ofstream f_out;
                utilities::GenFileStream(f_out, fileName, format, mpi);
                if(!mpi.m_exitFlag)
                {
                    if(io::_BINARY_ == format)
                    {
                        iSize_t strSize = m_fileHeader.size();
                        f_out.write((char*)&(strSize), sizeof(iSize_t));
                        f_out.write(m_fileHeader.c_str(), strSize);
                        f_out.write((char*)&nbrLabels, sizeof(iSize_t));
                        f_out.write((char*)&m_kMax, sizeof(iSize_t));
                        f_out.write((char*)m_vTable.data(), m_vTable.size()*sizeof(T));
                    }
                    else if(io::_TEXT_ == format)
                    {
                        f_out << m_fileHeader << "\n";
                        f_out << nbrLabels << "\n";
                        f_out << m_kMax << "\n";
                        for(kState_t k1=0; k1<m_kMax; ++k1)
                        {
                            std::string streamData = utilities::ToStream(m_vTable[k1]);
                            f_out << k1 << "\t" << k1 << "\t" << streamData <<"\n";
                        }
                    }
                    f_out.close();
                }
            }
            MPI_Barrier(mpi.m_comm);
            mpi.ExitFlagTest();
        }
        
        //!
        //! Read quadratic terms from a file
        //!
        void FromFile(
            const std::string fileName,     //!<    Name of file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class 
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                std::ifstream f_in;
                utilities::GenFileStream(f_in, fileName, format, mpi);
                if(!mpi.m_exitFlag)
                {
                    iSize_t nbrLabels = 0;
                    std::string fileHeader;
                    if(io::_BINARY_ == format)
                    {
                        iSize_t strSize;
                        f_in.read(reinterpret_cast<char*>(&strSize), sizeof(iSize_t));
                        fileHeader.resize(strSize);
                        f_in.read(&fileHeader[0], strSize);
                        if(m_fileHeader != fileHeader)
                        {
                            std::cerr << "\n\tFile read in error - wrong file header " << fileHeader << std::endl;
                            mpi.m_exitFlag = true;
                            f_in.close();
                            goto escape;
                        }
                        f_in.read(reinterpret_cast<char*>(&nbrLabels), sizeof(iSize_t));
                        f_in.read(reinterpret_cast<char*>(&m_kMax), sizeof(iSize_t));
                        this->Initialize(m_kMax, mpi);
                        f_in.read(reinterpret_cast<char*>(m_vTable.data()), m_vTable.size()*sizeof(T));
                    }
                    else if(io::_TEXT_ == format)
                    {
                        getline(f_in, fileHeader);
                        if(m_fileHeader != fileHeader)
                        {
                            std::cerr << "\n\tFile read in error - wrong file header " << fileHeader << std::endl;
                            mpi.m_exitFlag = true;
                            f_in.close();
                            goto escape;
                        }
                        f_in >> nbrLabels;
                        f_in >> m_kMax;
                        this->Initialize(m_kMax, mpi);
                        kState_t k1;
                        kState_t k2;
                        for(iSize_t i=0; i<m_vTable.size(); ++i)
                        {
                            f_in >> k1 >> k2;
                            utilities::FromStream(f_in, m_vTable[i]);
                        }
                    }
                    f_in.close();
                }
            }
            escape:
            {}
            mpi.ExitFlagTest();
            mpi.Sync(&m_kMax, 1, 0);
            if(0 != mpi.m_id)	// FOR ALL OTHER NODES
            {
                this->Initialize(m_kMax, mpi);
            }
            MPI_Barrier(mpi.m_comm);
            this->MpiSynchronize(0, mpi);
        }
    };
}   //  End namespace diagonalization
#endif
