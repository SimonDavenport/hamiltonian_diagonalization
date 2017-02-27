////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a wrapper for the FFTW library
//!
//!                    Copyright (C) Simon C Davenport
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

#ifndef _FOURIER_TRANSFORM_HPP_INCLUDED_
#define _FOURIER_TRANSFORM_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <fftw3.h>
#include "../general/dcmplx_type_def.hpp"
#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{
    void DiscreteFourierTransform1D(int N, int M, dcmplx* coefficients, dcmplx* outputs, int sign);
};  //  End namespace utilities
#endif
