////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		Header file for binary number tools
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
 
#ifndef _BINARY_NUMBER_TOOLS_HPP_INCLUDED_
#define _BINARY_NUMBER_TOOLS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "binomials.hpp"
#include <bitset>
#include <cstdint>
#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{
    namespace binary
    {
        ///////     STATIC CONST DECLARATIONS     //////////////////////////////////////
        static const uint64_t  m1  = 0x5555555555555555; //!<	Binary: 0101...
        static const uint64_t  m2  = 0x3333333333333333; //!<	Binary: 00110011..
        static const uint64_t  m4  = 0x0f0f0f0f0f0f0f0f; //!<	Binary:  4 zeros, 4 ones ...
        static const uint64_t  h01 = 0x0101010101010101; //!<	The sum of 256 to the power of 0,1,2,3...
        static const uint64_t  number1 = 1;			     //!<	The number 1 as an uint64_t

        //!
        //! A table of numbers in the B(2,6) De Bruijn squence, required for 
        //! binary calculation of logarithms (in base 2)
        //!
        static const int deBruijnBitSequence[64] = 
        {
            0,1,2,53,3,7,54,27,4,38,41,8,34,55,48,28,
            62,5,39,46,44,42,22,9,24,35,59,56,49,18,29,11,
            63,52,6,26,37,40,33,47,61,45,43,21,23,58,17,10,
            51,25,36,32,60,20,57,16,50,31,19,15,30,14,13,12,
        };
        
        static const uint64_t deBruijnMultiply = 0x022FDD63CC95386D;
        //!<    Required to address the B(2,6) De Bruijn squence
        static const uint64_t deBruijnShift = 58;
        //!<    Required to address the B(2,6) De Bruijn squence
        
        //////      FUNCTION DECLARATIONS       ////////////////////////////////////////
        int HammingWeight64(uint64_t x);
        int HammingWeight64Iterative(uint64_t x);
        uint64_t FirstBinaryHammingNumber64(const unsigned int hammingWeight);
        uint64_t IsolateLeftMostSetBit(uint64_t x);
        unsigned int Base2Log(uint64_t x);
        uint64_t NextHammingNumber64(const uint64_t x);
        uint64_t IndexHammingNumber64(const unsigned int hammingWeight,const uint64_t x);
        uint64_t GenerateHammingNumber(const unsigned int hammingWeight,const uint64_t index);   

	    ////////////////////////////////////////////////////////////////////////////////
        //! \brief Round down to the nearest power of 2 or in other words isolate the
        //! left-most bit
        //!
        //! \return Input rounded down to the nearest power of 2
        ////////////////////////////////////////////////////////////////////////////////
	    template <int N>
	    std::bitset<N> IsolateLeftMostSetBit(
	        std::bitset<N>  x)     //!<    A number to find the left-most set bit of
	    {
            int shift = 1;
            while(shift<=N/2)
            {
                x |= x >> shift;
                shift*=2;
            }
            ++x;
            return x >> number1;
        }
    }   //  End namespace binary
}   // End namespace utilities
#endif
