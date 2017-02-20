////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 27/03/2015              
//!                                                                             
//!	 \file
//!     This file generates a set of two-body pseudopotential coefficients
//!     for the Coulomb interaction, implemented in the sphere geometry        
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

#ifndef _COULOMB_PSEUDOPOTENTIALS_HPP_INCLUDED_
#define _COULOMB_PSEUDOPOTENTIALS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../../utilities/general/orbital_and_state_defs.hpp"
                                            //  Declaration of iSize_t
#include "../../utilities/mathematics/binomial_table.hpp"   
                                            //  For calculation of Binomials
#include <vector>                           //  For std::vector
#include <iostream>                         //  For std::cerr

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    std::vector<double> GenCoulombPseudopotentials(
        const iSize_t nbrOrbitals,
        const iSize_t llIndex,
        const double mulFactor);
    
}   //  End namespace diagonalization 

#endif
