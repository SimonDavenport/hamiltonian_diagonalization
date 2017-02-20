////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 06/11/2014
//!
//!  \file
//!		This file defines some types that are generically used to represent
//!     fermion orbital labales and binary representations of a fermion
//!     Fock space
//!
//!                    Copyright (C) 2014 Simon C Davenport
//!
//!		This program is free software: you can redistribute it and/or modify
//!		it under the terms of the GNU General Public License as published by
//!		the Free Software Foundation, either version 3 of the License,
//!		or (at your option) any later version.
//!
//!		This program is distributed in the hope that it will be useful, but
//!		WITHOUT ANY WARRANTY; without even the implied warranty of 
//!		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
//!		General Public License for more details.
//!
//!		You should have received a copy of the GNU General Public License
//!		along with this program. If not, see <http://www.gnu.org/licenses/>.
//! 
//////////////////////////////////////////////////////////////////////////////// 

#ifndef _ORBITAL_AND_STATE_DEFS_HPP_INCLUDED_
#define _ORBITAL_AND_STATE_DEFS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <cstdint>  //  for "uint" types
#include <limits>   //  for std::numeric_limits

namespace diagonalization
{

///////     VARIABLE TYPE DECLARATIONS     /////////////////////////////////////

typedef uint64_t fock_t;
//!<    Define a Fock occupation basis for fermions, simply a bitwise encoding
//!     of the set of orbital occupations (up to 64 orbitals can be addressed)   

typedef uint16_t kState_t;
//!<    Define an unsigned quantum number index type (a 16 bit unsigned int)

typedef int16_t mState_t;
//!<    Define a signed quantum number index type (a 16 bit signed int)

typedef uint32_t iSize_t;
//!<    Define a size type to store parameters such as particle number and 
//!     total orbital number, and general unsigned ints (a 32 bit int)

////////////////////////////////////////////////////////////////////////////////
//! \brief Set the type of index to address the sparse matrix with (e.g. uint_64_t)
//!
//! NOTE: using int rather than long int is faster and saves memory space, 
//! but limits the maximum matrix size to 2^31 -1 i.e. 2,147,483,647
//!
//! NOTE 2: When used during FORTRAN matrix-vector routine, only int type
//! can be passed.
//!
////////////////////////////////////////////////////////////////////////////////
typedef int crsIndex_t;

////////////////////////////////////////////////////////////////////////////////
//! \brief Error code to return if a Fock state is found beyond a certian limit
//!
////////////////////////////////////////////////////////////////////////////////

static const fock_t _SEARCH_ERROR_CODE_ = std::numeric_limits<fock_t>::max();

}   //  End namespace diagonalization 

#endif

