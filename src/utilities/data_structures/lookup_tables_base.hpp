////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                           
//!                                                                             
//!                      \date Last Modified: 02/02/2015                        
//!                                                                             
//!	 \file
//!     This file contains the base class implementation for lookup tables
//!                                                  
//!                    Copyright (C) 2015 Simon C Davenport
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

#ifndef _LOOKUP_TABLES_BASE_HPP_INCLUDED_
#define _LOOKUP_TABLES_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../general/orbital_and_state_defs.hpp"
#include "../wrappers/mpi_wrapper.hpp"
#include <vector>

#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A class to contain look-up table values and their associated
    //! momentum conservation laws
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    class LookUpTables
    {
        protected:
        
        std::vector<T> m_vTable;            //!<    Look-up table data 
        
        std::vector<kState_t> m_kTable;     //!<    Store a look-up table of k1 values 
                                            //!     satisfying a conservation law

        kState_t m_kMax;                    //!<    Maximum k-value
        
        //  Function to return the dimension of the v table
        virtual iSize_t CalculateDim(const kState_t kMax) const=0;   
        
        public:
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Default constructor
        //!
        ////////////////////////////////////////////////////////////////////////////////

        LookUpTables()
        :
            m_kMax(0)
        {}
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Destructor
        //!
        ////////////////////////////////////////////////////////////////////////////////

        ~LookUpTables()
        {
            this->Clear();
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Initialize memory allocations
        //!
        ////////////////////////////////////////////////////////////////////////////////

        void Initialize(
            const kState_t kMax,                //!<    No. k states
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class 
        {
            m_kMax = kMax;

            iSize_t dim = this->CalculateDim(kMax);
            
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            { 
                utilities::cout.AdditionalInfo()<<"\n\t- INITIALIZING LOOKUP TABLE "<<std::endl;
            
                utilities::cout.DebuggingInfo()<<"\n\t- ALLOCATING "<<2*sizeof(T)*dim/(1024.0*1024.0)
                         <<" mb to store look-up table of "<<dim<<" lookup terms."<<std::endl;
            }
            
            m_kTable.resize(dim);
            m_vTable.resize(dim);
        }
       
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Reset the memory allocations.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        void Clear()
        {
            m_kTable.clear();
            m_vTable.clear();
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //  Get methods: default implementation here, optionally  
        //  implemented in derived classes
        
        virtual void GetK1(kState_t* kRetrieveBuffer,iSize_t& nbrK1,const kState_t k2,const kState_t k3,const kState_t k4) const{};

        virtual T GetVkkkk(const kState_t k1,const kState_t k2,const kState_t k3,const kState_t k4) const{ return 0.0;};
        
        virtual void GetK1(kState_t* kRetrieveBuffer,iSize_t& nbrK1,const kState_t k2) const{};

        virtual T GetEkk(const kState_t k1,const kState_t k2) const{ return 0.0;};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Only allow for a single k value to be returned
        //!
        ////////////////////////////////////////////////////////////////////////////////

        iSize_t GetMaxKCount() const
        {
            return 1;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Get the dimension of the quartic term array
        //!
        ////////////////////////////////////////////////////////////////////////////////

        iSize_t GetDimension() const
        {
            return m_vTable.size();
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Return the address of the class momentum table so that 
        //! it can be set externally.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        std::vector<kState_t>* GetKTable()
        {
            return &m_kTable;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Return the address of the class quadratic term array so 
        //! that it can be set externally.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        std::vector<T>* GetVTable()
        {
            return &m_vTable;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Given an overall coefficient g multiplying all of the 
        //! quartic terms, this function incorporates the coefficient into 
        //! all table values.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        void SetCoefficient(
            const T coefficient)   //!<    Coefficient value
        {
            for(auto& it : m_vTable)
            {
                it *= coefficient;
            }
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Synchronize table contents with a given node
        //!
        ////////////////////////////////////////////////////////////////////////////////

        void MpiSynchronize(
            const iSize_t syncId,               //!<    Node to synchronize with
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            mpi.Sync(&m_vTable,syncId);
            mpi.Sync(&m_kTable,syncId);
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Write the table to a binary file.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        void TableToFile(
            const std::string fileName, //!<    Name of file
            utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
            const
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                std::ofstream f_out;
            
                f_out.open(fileName.c_str(), std::ios::binary);

                if(f_out.is_open())
                {
                    iSize_t dim = m_vTable.size();
                
                    f_out.write((char*)&dim,sizeof(iSize_t));
                    f_out.write((char*)m_vTable.data(),dim*sizeof(T));
                    f_out.close();
                    
                    utilities::cout.AdditionalInfo()<<"DONE (see e.g. "<<fileName<<")"<<std::endl;
                }
                else
                {
                    std::cerr<<"ERROR in TableToFile(). File "<<fileName<<" not found!"<<std::endl;
                    mpi.m_exitFlag=true;
                }
            }
            
            //	Barrier waits for all nodes to catch up to this point
            MPI_Barrier(mpi.m_comm);
            
            mpi.ExitFlagTest();
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Read the table in from an existing binary file.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        void TableFromFile(
            const std::string fileName, //!<    Name of file
            utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                std::ifstream f_in;
                
                f_in.open(fileName.c_str(), std::ios::binary);

                if(f_in.is_open())
                {
                    iSize_t dim = 0;
                
                    f_in.read(reinterpret_cast<char*>(&dim),sizeof(iSize_t));
                    
                    //  Allocate memory to store the table

                    m_vTable.resize(dim);
                    
                    f_in.read(reinterpret_cast<char*>(m_vTable.data()),dim*sizeof(dcmplx));
                    f_in.close();
                    utilities::cout.SecondaryOutput()<<"DONE ";
                }
                else
                {
                    std::cerr<<"ERROR in TableFromFile(). File "<<fileName<<" not found!"<<std::endl;
                    mpi.m_exitFlag=true;
                }        
            }
            
            mpi.ExitFlagTest();

            //  MPI sync the tables

            mpi.Sync(&m_vTable,0);
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End diagonalization namespace

#endif

