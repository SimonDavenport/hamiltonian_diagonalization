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

#ifndef _TERM_HASH_TABLES_BASE_HPP_INCLUDED_
#define _TERM_HASH_TABLES_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../general/orbital_and_state_defs.hpp"
#include "../data_structures/multi_key_hash.hpp"
#include "../wrappers/mpi_wrapper.hpp"
#include <vector>

namespace diagonalization
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A class to contain term table values of the form 
    //! V_{k1,k2,k3,k4}. This version uses an underlying hash table
    //! data structure.
    ////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    class TermHashTables
    {
        protected:
        utilities::MultiHashMap<T> m_vTable;    //!< Look-up table data
        utilities::MultiHashMultiMap<kState_t> m_kTable;
                                            //!< Store a look-up table of k1 values 
                                            //!  satisfying a conservation law   
        kState_t m_kMax;                    //!<    Maximum index-value
        public:

        //!
        //! Default constructor
        //!
        TermHashTables()
        :
            m_kMax(0)
        {}

        //!
        //! Destructor
        //!
        ~TermHashTables(){}

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

        virtual void GetK1(kState_t* kRetrieveBuffer, iSize_t& nbrK1, const kState_t k2,
                           const kState_t k3, const kState_t k4) const{};
        virtual T GetVkkkk(const kState_t k1, const kState_t k2, const kState_t k3,
                           const kState_t k4) const{return 0.0;};
        virtual void GetK1(kState_t* kRetrieveBuffer, iSize_t& nbrK1, const kState_t k2) const{};
        virtual T GetEkk(const kState_t k1, const kState_t k2) const{return 0.0;};
        virtual iSize_t GetMaxKCount() const{return 1;};
        virtual void ToFile(const std::string fileName, const std::string format, 
                                 utilities::MpiWrapper& mpi)=0;
        virtual void FromFile(const std::string fileName, const std::string format,
                                   utilities::MpiWrapper& mpi)=0;
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
        //! Output the table to a file
        //!
        void ToFileBase(
            const std::string fileName,         //!<    File name
            const std::string format,           //!<    Format of file (e.g. "binary", "text")
            const iSize_t nbrLabels,            //!<    Number of quantum number labels
            utilities::MpiWrapper& mpi)         //!<    Instance of the MPI wrapper class
        {
            m_vTable.ToFile(fileName, format, nbrLabels, mpi);
        }

        //!
        //! Import the term table from a file
        //!
        void FromFileBase(
            const std::string fileName,         //!<    File name
            const std::string format,           //!<    Format of file (e.g. "binary", "text")
            utilities::MpiWrapper& mpi)         //!<    Instance of the MPI wrapper class
        {
            m_vTable.FromFile(fileName, format, mpi);
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
        void SetQuarticFromArray(
            C* quarticArray)    //!<    Class container for quartic term array
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
                        m_vTable.Insert(utilities::Key(k1, k2, k3, k4)) = quarticArray->GetVkkkk(k1, k2, k3, k4);
                    }
                }
            }
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
                m_vTable.Insert(utilities::Key(k2, k2)) = quadraticArray->GetEkk(k1, k2);
            }
        }
    };
}   //  End namespace diagonalization 
#endif
