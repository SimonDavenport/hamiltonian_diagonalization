////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 18/09/2014
//!
//!  \file
//!		This file contains a function template for a modulo operator that deals
//!     properly with negative numbers. 
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

#ifndef _MODULO_HPP_INCLUDED_
#define _MODULO_HPP_INCLUDED_

#include <cmath> //  For definition of std::abs

namespace utilities
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A template function that returns the result a modulo b, for
    //! both positive and negative values of a. 
    //!
    //! If a is negative and b is positive then e.g. -5 % 4 gives -1 because
    //! it acts on the numerical part of the binary, not the sign:
    //!
    //! -5 is 1|101         and    4 is 0|100
    //!       ^ sign part
    //!
    //! then 1|101 % 0|100 becomes a % operation on the part not representing 
    //! the sign, so we get -(5%4), which gives -1. 
    //!
    //! We really wanted it to give 3, which is obtained by calculating
    //! 4 - abs(-1) % 4 = 3
    //!
    //! NOTE: negative values of b are not treated correctly here. 
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    template <typename T>
    T Modulo(
        const T a,    //!<    First argument in modulo function
        const T b)    //!<    Second argument in modulo function
    {
        return a >= 0 ? a % b : ( b - (T)std::abs ( a%b ) ) % b;
    }

}   //  End namespace utilities

#endif
