////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                  
//!                                                                             
//!	 \file
//!     This file defines a class to store data structures describing 
//!     coefficients in a Hamiltonian
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

#ifndef _LOOKUP_TABLES_HPP_INCLUDED_
#define _LOOKUP_TABLES_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/general/dcmplx_type_def.hpp"
#include "../../utilities/data_structures/lookup_tables_base.hpp"
#include "../../utilities/data_structures/lookup_hash_tables_base.hpp"
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif    

namespace diagonalization
{
    class QuadraticLookUpTables : public LookUpTables<dcmplx>
    {
        private:
        iSize_t CalculateDim(const kState_t kMax) const;
        public:
        void GetK1(kState_t* kRetrieveBuffer,iSize_t& nbrK1,const kState_t k2) const;
        dcmplx GetEkk(const kState_t k1,const kState_t k2) const;
    };
    
    class QuarticLookUpTables : public LookUpTables<dcmplx>
    {
        private:
        iSize_t CalculateDim(const kState_t kMax) const;
        public:
        void GetK1(kState_t* kRetrieveBuffer,iSize_t& nbrK1,const kState_t k2,const kState_t k3,const kState_t k4) const;
        dcmplx GetVkkkk(const kState_t k1,const kState_t k2,const kState_t k3,const kState_t k4) const;
    };
    
    class QuadraticLookUpHashTables : public LookUpHashTables<dcmplx>
    {
        public:
        iSize_t GetMaxKCount() const;
        void GetK1(kState_t* kRetrieveBuffer,iSize_t& nbrK1,const kState_t k2) const;
        dcmplx GetEkk(const kState_t k1,const kState_t k2) const;
    };
    
    class QuarticLookUpHashTables : public LookUpHashTables<dcmplx>
    {
        public:
        iSize_t GetMaxKCount() const;
        void GetK1(kState_t* kRetrieveBuffer,iSize_t& nbrK1,const kState_t k2,const kState_t k3,const kState_t k4) const;
        dcmplx GetVkkkk(const kState_t k1,const kState_t k2,const kState_t k3,const kState_t k4) const;
    };
}   //  End diagonalization namespace
#endif