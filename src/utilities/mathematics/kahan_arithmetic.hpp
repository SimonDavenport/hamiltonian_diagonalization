////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This header file implements a  class template that implements the
//!     Kahan summation algorithm on all values added to the object with += 
//!
//!     The Kahan Summation algorithm is implemented in order to reduce
//!     numerical errors.
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

#ifndef _KAHAN_ARITHMETIC_HPP_INCLUDED_
#define _KAHAN_ARITHMETIC_HPP_INCLUDED_

#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{
    //!
    //! A class template to implement the Kahan summation algorithm
    //!
    template<typename T>
    struct KahanAccumulation
    {
        T m_sum;                //!<    Current accumulated sum 
        T m_correction;         //!<    Correction term
        //!
        //! Default constructor initialises sum to zero
        //! 
        KahanAccumulation()
        :
            m_sum(0),
            m_correction(0)
        {}
    };

    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief Defines a global  template += operator overload for the Kahan 
    //! accumulator class. Call as simply a += b where a is an existing 
    //! KahanAccumulation object and b is the numerical value to be added to 
    //! the sum (should be of the same type as the class type)
    //!
    //! \return address to the updated KahanAccumulation object
    //!
    //////////////////////////////////////////////////////////////////////////////// 
    template<typename T>
    KahanAccumulation<T>& operator+=(KahanAccumulation<T>& a, const T& b)
    {
        T tempVal = b - a.m_correction;
        T tempSum = a.m_sum + tempVal;
        a.m_correction = (tempSum - a.m_sum) - tempVal;
        a.m_sum = tempSum;
        return a;
    }
}   //  End namespace utilities
#endif
