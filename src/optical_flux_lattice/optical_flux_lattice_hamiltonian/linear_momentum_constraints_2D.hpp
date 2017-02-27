////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                           
//!                                                                             
//!	 \file
//!     This file defines functions that impose 2D total linear momentum
//!     constraints on a given fock state
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

#ifndef _LINEAR_MOMENTUM_CONSTRAINTS_HPP_INCLUDED_
#define _LINEAR_MOMENTUM_CONSTRAINTS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/mathematics/modulo.hpp"
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../utilities/mathematics/binary_number_tools.hpp"

namespace diagonalization
{
    bool TestLinearMomentumSector2(const fock_t state, const kState_t kxTot,
                                   const kState_t kyTot, const kState_t dimX, 
                                   const kState_t dimY); 
    bool TestLinearMomentumSector1(const fock_t state, const kState_t kyTot,
                                   const kState_t dimY);
    kState_t StateAddition(const kState_t k1, const kState_t k2,
                           const kState_t dimX1, const kState_t dimY1,
                           const kState_t dimX2, const kState_t dimY2);
    kState_t StateSubtraction(const kState_t k1, const kState_t k2,
                              const kState_t dimX1, const kState_t dimY1,
                              const kState_t dimX2, const kState_t dimY2);
    kState_t R60(const kState_t k, const int r, const kState_t dimX,
                 const kState_t dimY);
    kState_t InverseR60(const kState_t k, const int r, const kState_t dimX,
                        const kState_t dimY);
}   //  End namespace diagonalization
#endif
