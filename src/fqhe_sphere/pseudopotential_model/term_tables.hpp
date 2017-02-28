////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store data structures describing 
//!     terms in a Hamiltonian
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

#ifndef _TERM_TABLES_HPP_INCLUDED_
#define _TERM_TABLES_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/data_structures/term_tables_base.hpp"
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    class QuadraticTermTables : public TermTables<double>
    {
        private:
        iSize_t CalculateDim(const kState_t kMax) const;
        public:
        void GetK1(kState_t* kRetrieveBuffer, iSize_t& nbrK1, const kState_t k2) const;
        double GetEkk(const kState_t k1, const kState_t k2) const;
        void ToFile(const std::string fileName, const std::string format, 
                    utilities::MpiWrapper& mpi) const;
        void FromFile(const std::string fileName, const std::string format,
                      utilities::MpiWrapper& mpi);
    };
    
    class QuarticTermTables : public TermTables<double>
    {
        private:
        iSize_t CalculateDim(const kState_t kMax) const;
        public:
        void GetK1(kState_t* kRetrieveBuffer, iSize_t& nbrK1, const kState_t k2, 
                   const kState_t k3, const kState_t k4) const;
        void SetK1(const kState_t k1, const kState_t k2, 
                   const kState_t k3, const kState_t k4);
        double GetVkkkk(const kState_t k1, const kState_t k2, const kState_t k3, 
                        const kState_t k4) const;
        void SetVkkkk(const double Vkkkk, const kState_t k1, const kState_t k2, 
                      const kState_t k3, const kState_t k4);
        void ToFile(const std::string fileName, const std::string format, 
                    utilities::MpiWrapper& mpi) const;
        void FromFile(const std::string fileName, const std::string format,
                      utilities::MpiWrapper& mpi);
    };
}   //  End namespace diagonalization
#endif
