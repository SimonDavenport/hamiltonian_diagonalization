////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 05/04/2014
//!
//!  \file
//!		This file contains some efficient algorithms for evaluating the Hamming
//!     weight of a binary number - i.e. the number of bits that are set to 1.
//!     Also, there are algorithms to generate a list of binary numbers with
//!		the same Hamming weight
//!
//!                    Copyright (C) 2014 Simon C Davenport
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

#include "binary_number_tools.hpp" 

namespace utilities
{

namespace binary
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
	
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Hamming weight algorithm with minimal number of arithmetic operations
    //!
    //! THIS ALGORITHM IS ADAPTED FROM THE WIKIPEDIA ARTICLE: 
    //!     http://en.wikipedia.org/wiki/Hamming_weight
    //!
    //! x86 implementation requires 24 instructions
    //!
    //! Example:        
    //!
    //! x = 00 11 01 11 11 01 11 00 10 11 11 11 11 01 01 01
    //!
    //! (x>>1) & 010101010101... = 00 01 00 01 01 00 01 00 01 01 01 01 01 00 00 00
    //!
    //! This function maps bit pairs in the following way:
    //!
    //! 11 -> 01  10  -> 01  01 -> 00  00 -> 00 
    //!
    //! Then x -= (x>>1) & 010101010101... gives
    //!
    //!      00 10 01 10 10 01 10 00 01 10 10 10 10 01 01 01
    //! i.e. 0  2  1 2  2  1  2  0  1  2  2  2  2  1  1  1
    //!
    //! This gives the count of sets bits from each original bit pair
    //!
    //! (x & m2) + ((x >> 2) & m2) simply adds the numbers stored in pairs 
    //! of 2-bits using that m2 = 001100110011...)
    //! 
    //! NOTE (x + (x >> 2)) & m2 won't work because there can be 4 set bits
    //! in each 4-bit section - this is represetned by 0100, but this would
    //! by masked off by m2
    //!     
    //!      00 10 00 10 00 01 00 00 00 10 00 10 00 01 00 01
    //!    + 00 00 00 01 00 10 00 10 00 01 00 10 00 10 00 01
    //!      0010  0011  0011  0010  0011  0100  0011  0010
    //! i.e.    2     3     3     2     3     4     3     2
    //!
    //! The remaining steps are in the same spirit
    //!
    //!  (x + (x >> 4)) & 0000111100001111...
    //!
    //!      0010 0011 0011 0010 0011 0100 0011 0010
    //!    + xxxx 0010 0011 0011 0010 0011 0100 0011
    //!      xxxx 0101 0110 0101 0101 0111 0111 0101
    //!
    //! Then take & 0000111100001111...
    //!
    //!       00000101 00000101 00000111 00000101
    //! i.e.  2+3=5    2+3=5    3+4=7    2+3=5     
    //!
    //! Next we use a trick of multiplying by h01 = 0x0101010101010101
    //!
    //! x * h01 = 0001011000010110000101100001011000010110000100010000110000000101
    //!
    //! Finally, right shift by 56 (i.e. count only the left 8 bits)
    //!
    //! return 00010110 = 22=5+5+7+5
    //!
    //! \return Hamming weight (number of 1 bits in the number)
    //!
    ////////////////////////////////////////////////////////////////////////////////

    int HammingWeight64(
        uint64_t x)        //!<    Number to find the Hamming weight of
    {
        if(0==x) return 0;

        x -= (x >> 1) & m1;             //  Put count of each 2 bits into those 2 bits
        x = (x & m2) + ((x >> 2) & m2); //  Put count of each 4 bits into those 4 bits 
        x = (x + (x >> 4)) & m4;        //  Put count of each 8 bits into those 8 bits 
  
        return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
	
	////////////////////////////////////////////////////////////////////////////////
    //! \brief Iterative Hamming weight algorithm.
    //!
    //! This algorithm is more efficient when the number of bits counted 
	//!	is known to be very small.
    //!
    //! x86 implementation requires 13 instructions
    //!
    //! Example:
    //!
    //! x = 00 00 01 10 00 00 00 10 00 00 01 00 00 00 01 00
    //!
    //! x-1 = 00 00 01 10 00 00 00 10 00 00 01 00 00 00 00 11
    //!
    //! x &= x-1 gives 00 00 01 10 00 00 00 10 00 00 01 00 00 00 00 00
    //!
    //! x-1 = 00 00 01 10 00 00 00 10 00 00 01 11 11 11 11 11
    //!
    //! x &= x-1 gives 00 00 01 10 00 00 00 10 00 00 00 00 00 00 00 00
    //!
    //! x-1 = 00 00 01 10 00 00 00 01 11 11 11 11 11 11 11 11
    //!
    //! x &= x-1 gives 00 00 01 10 00 00 00 00 00 00 00 00 00 00 00 00
    //!
    //! etc.
    //!
    //! Notice that the x &= x-1 masks off all bits to the right of the right
    //! -most set bit. This operation needs to be completed M times, where M
    //! is the number of set bits
    //!
    //! \return Hamming weight (number of 1 bits in the number)
    //!
    ////////////////////////////////////////////////////////////////////////////////
 
    int HammingWeight64Iterative(
        uint64_t x)        //!<    Number to find the Hamming weight of
    {
        if(0==x) return 0;
    
        unsigned int i=1;
		
		while (x &= (x-binary::number1)) ++i;
		
        return i;  
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Round down to the nearest power of 2 or in other words isolate the
    //! left-most bit
    //!
    //! x86 implementation requires 22 instructions
    //!
    //! Example:
    //!
    //! x = 00 00 01 10 00 00 00 10 00 00 01 00 00 00 01 00
    //!
    //! x |= x >> 1 gives 00 00 01 11 00 00 00 11 00 00 01 10 00 00 01 10
    //!
    //! x |= x >> 2 gives 00 00 01 11 11 00 00 11 11 00 01 11 10 00 01 11
    //!
    //! x |= x >> 4 gives 00 00 01 11 11 11 11 11 11 11 11 11 11 11 11 11
    //! 
    //! The remaining operations do not change this arrangement. Then the
    //! last steps are to add on and right shift by 1, thereby isolating the
    //! left-most bit    
    //!
    //! \return Input rounded down to the nearest power of 2
    //!
    ////////////////////////////////////////////////////////////////////////////////
	
	uint64_t IsolateLeftMostSetBit(
	    uint64_t  x)     //!<    A number to find the left-most set bit of
	{
        x |= x >> 1; 
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        x |= x >> 32;
        ++x;

        return x >> binary::number1;
    }
 
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the base-2 log of a number (assuming that the input is 
    //! a power of 2). Alternatively, determines the position of a single set bit
    //! counting from the right
    //!
    //! x86 implementation requires 5 instructions
    //!
    //! \return The log_2 of the input
    //!
    ////////////////////////////////////////////////////////////////////////////////
        
    unsigned int Base2Log(
	    uint64_t  x)     //!<    A number to find log of that is a power of 2
    {
        //  deBruijnBitSequence has a static definition in binary_number_tools.h
        return deBruijnBitSequence[(uint64_t)(x * (uint64_t)deBruijnMultiply) >> deBruijnShift];
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate the "first" binary number of with a given Hamming weight
    //!
    //! First means that it's of the form ...0001...1111
    //!
    //! x86 implementation requires 5 instructions
    //!
    //! \return The uint64_t corresponding to the binary ...0001...1111
    //!
    ////////////////////////////////////////////////////////////////////////////////
	
	uint64_t FirstBinaryHammingNumber64(
	    const unsigned int hammingWeight)     //!<    Hamming weight of the number (i.e no. 1s in it)
	{
		return ((binary::number1 << hammingWeight) - binary::number1);
	}
	
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate the next lexicographic permutation of numbers with the
	//!	same Hamming weight
	//!
	//!	Algorithm derived from:
	//!
	//!	http://www.geeksforgeeks.org/next-higher-number-with-same-number-of-set-bits/
	//!
	//!	EXAMPLE:
	//!
	//!	10011100
	//!	00011100 - right most string of 1's in x
	//!	00000011 - right shifted pattern except left most bit ------> [A]
	//!	00010000 - isolated left most bit of right most 1's pattern
	//!	00100000 - shiftleft-ed the isolated bit by one position ------> [B]
	//!	10000000 - left part of x, excluding right most 1's pattern ------> [C]
	//!	10100000 - add B and C (OR operation) ------> [D]
	//!	10100011 - add A and D to get the final result
    //!
    //! x86 implementation requires 15 instructions
    //!
    //! Example:
    //!
    //! x = 00 00 00 00 00 00 00 00 00 00 00 00 00 10 01 11 00
    //!
    //! -x gives the 2s complement (flip bits, add 1)
    //!
    //! temp = x & -x = 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 01 00
    //!                 
    //! x + temp = 00 00 00 00 00 00 00 00 00 00 00 00 00 10 10 00 00  [D]
    //!              
    //! x ^ (x + temp) = 00 00 00 00 00 00 00 00 00 00 00 00 00 00 11 11 00
    //!
    //! x ^ (x + temp) /temp = 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 11 11
    //!
    //! previous line >> 2 = 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 11
    //!
    //! OR x ^ (x + temp) >> 2+log_2 temp gives the same result
    //!
    //! final output is 00 00 00 00 00 00 00 00 00 00 00 00 00 10 10 00 11
    //!
    //! \return The lexicographically next binary number of the same Hamming weight
    //!
    ////////////////////////////////////////////////////////////////////////////////

	uint64_t NextHammingNumber64(
	    const uint64_t x)	//!<    A number with some bitwise representation
	{   
	    if(0==x)    return x;
	
		uint64_t temp = (x & -x);  //  Isolates the right-most bit
		
		//return ((x ^ (x + temp)) >> (2+utilities::binary::Base2Log(temp))) | (x + temp);
		return ((x ^ (x + temp)) / temp) >> 2 | (x + temp);
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate the index in the list of lexicographically ordered Hamming
    //! numbers assuming all possible numbers are included in order
    //!
    //! \return The position of the input in a lexicographically ordered list of 
    //! hamming numbers
    //!
    ////////////////////////////////////////////////////////////////////////////////

	uint64_t IndexHammingNumber64(
	    const unsigned int hammingWeight,     //!<    Hamming weight of the set
	    const uint64_t x)	                  //!<    A number that we want to find the index of
	{
		uint64_t output = 0;
		uint64_t temp = x;
		
		for(int p = hammingWeight; p>0; --p)
		{
    		uint64_t leftMostBit = IsolateLeftMostSetBit(temp);         
    		//  Isolates the left-most bit

            //  Find the index of the left-most set bit
		    int m = Base2Log(leftMostBit);
		    
		    //std::cout<<"\t\ttemp="<<temp<<" IsolateLeftMostSetBit(temp) = "<<leftMostBit<<" power = "<<m<<std::endl;
		    //getchar();
		    
		    if(p<=m)
		    {
                //  Add binomial(m,p) from a look-up table with additional conditions

                //output += boost::math::binomial_coefficient<double>(m,p);
                
                output += utilities::BinomialFromTable(m,p);
                
                //if(boost::math::binomial_coefficient<double>(m,p) != utilities::BinomialFromTable(m,p))
                //{
                //    std::cerr<<" ERROR on node "<<utilities::mpi.m_id<<" with m = "<<m<<" p= "<<p<<std::endl;
                //}

                //std::cout<<"\t\tm="<<m<<" p="<<p<<" output = "<<output<<std::endl;

		    }
            
            //  Remove the current left-most bit
            
		    temp ^= leftMostBit;
		}

		return output;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate the Hamming number with the given hamming weight and at
    //! a given position in a lexicographic table (without having to generate 
    //! the full table of lexigographic permutations). This assumes a full list
    //! of numbers is going to be indexed
    //!
    //! \return The Hamming number at the index returnIndex
    //!
    ////////////////////////////////////////////////////////////////////////////////
	
    uint64_t GenerateHammingNumber(
        const unsigned int hammingWeight,   //!<    Hamming weight of the set
        const uint64_t index)               //!<    Index to return the Hamming number for
    {
        if(0 == index)
        {
            return FirstBinaryHammingNumber64(hammingWeight);
        }
        else
        {           
            //  Determine the exact binary Hamming number corresponding to this
            //  index (assuming lexicographic ordering)

            uint64_t returnVal = 0;
            uint64_t currIndex = 0;
            
            //std::cout<<"\t\t index = "<<index<<std::endl;
            
            for(int p = hammingWeight-1;p>=0;--p)
            {
                //  Determine the position of the left-most bit
                //  then next left most bit etc.
                
                int m = p;
                uint64_t incrementIndex;

                do
                {
                    //incrementIndex = boost::math::binomial_coefficient<double>(m,p);
                
                    incrementIndex = utilities::BinomialFromTable(m,p);
                
                    //if(boost::math::binomial_coefficient<double>(m,p) != utilities::BinomialFromTable(m,p))
                    //{
                    //    std::cerr<<" ERROR on node "<<utilities::mpi.m_id<<" with m = "<<m<<" p= "<<p<<std::endl;
                    //}
                    
                    currIndex += incrementIndex;
                    ++m;
                    
                    //std::cout<<"\t\tm="<<m<<" p="<<p<<" currIndex = "<<currIndex<<std::endl;
                    //getchar();
                }
                while(index >= currIndex);
                
                //  We've gone one step too far here, so go back one step
                
                currIndex -= incrementIndex;
                
                //std::cout<<"\t\tm="<<m<<" p="<<p<<" currIndex= "<<currIndex<<std::endl;
                
                //  Set the next left-most bit in the return value
                
                returnVal |= (binary::number1 << (m-1));
            }

            return returnVal;
        }
    }
  
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   // End namespace binary

}   // End namespace utilities
