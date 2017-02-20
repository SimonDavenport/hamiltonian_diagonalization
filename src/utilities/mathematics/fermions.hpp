////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 12/12/2014
//!
//!  \file
//!		Functions to perform binary encodings of fermion operator commutations
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

#ifndef _FERMIONS_HPP_INCLUDED_
#define _FERMIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <cstdint>                  //  For uint64_t
#include "binary_number_tools.hpp"  //  For Hamming weight algorithms

#if _DEBUG_
#include "../general/debug.hpp"
#include <bitset>
#endif

namespace utilities
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Apply an operator c^+_k1 c^+_k2 c_k3 c_k4 to a given state 
//! (with k1+k2=k3+k4)
//!
//! There will be only one possible resulting state since we are dealing with 
//! fermions (bosons is more complicated, since in that case there can be multiple
//! final states). 
//!
//!	NOTE: the states are assumed to be written in ascending order of orbital labels
//!
//! NOTE: this implementation assumes a 64-bit architecture
//!
//! \return The index of the final state.
//!
////////////////////////////////////////////////////////////////////////////////

inline uint64_t Evaluate_CdCdCC(
    int&  coefficient,              //!<    Address to store the sign of the output state  
    const uint64_t k1Occupied,      //!<    State with first orbital index occupied
    const uint64_t k2Occupied,      //!<    State with second orbital index occupied  
    const uint64_t k3Occupied,      //!<    State with third orbital index occupied   
    const uint64_t k4Occupied,      //!<    State with fourth orbital index occupied  
    const uint64_t state,           //!<    The state to be operated on
    const uint64_t highestOrbital)  //!<    Highest possible orbital
{
	//	Store the following results for subsequent use. These values store
	//	e.g. ...00000010000...
	//				  ^ - corresponds to the k1,k2,k3 or k4 orbital
	//	Note that we denote ...0001 as the 0th orbital

    #if _DEBUG_
        std::cout<<"k1: "<<std::bitset<64>(k1Occupied)<<" k2: "<<std::bitset<64>(k2Occupied)<<" k3: "<<std::bitset<64>(k3Occupied)<<" k4: "<<std::bitset<64>(k4Occupied)<<" state: "<<std::bitset<64>(state)<<std::endl;
    #endif

    //  Check conditions for which the result will be zero
    //  NOTE: lowest state is const static uint64_t 1
    if( 
    (k3Occupied>highestOrbital)        //  Orbitals higher than highestOrbital are   
    || (k4Occupied>highestOrbital)     //  not occupied                                               
    || ((state & k3Occupied)==0)         //  If the state does not already contain 
    || ((state & k4Occupied)==0)         //  both a k3 and k4 orbitals then c_k3 c_k4
                                         //  acts to give zero 
    || (k1Occupied & k2Occupied) ||  (k3Occupied & k4Occupied)     //  If a c or c^+ is repeated, then we get zero
    || (!(k1Occupied & k3Occupied) && !(k1Occupied & k4Occupied) && ((state & k1Occupied)!=0))   //	If c^+_k1 c^+_k2 appears it vanishes, 
	|| (!(k2Occupied & k3Occupied) && !(k2Occupied & k4Occupied) && ((state & k2Occupied)!=0))   // but only if k1!=k3 or k4 and k2!=k3 or k4. 
    )
    {   
        coefficient = 0;
        
        #if _DEBUG_
        std::cout<<"return 0"<<std::endl;
        #endif
        
        return 0;    //  Default return, but it won't be used   
    }

    //  If the input state is not zero then the output state is related to 
    //  it by a simple set of bitshift operations.
    //  We also need to keep track of the sign of the output state
	
	//	Declare some working variables
	
	uint64_t outState;
	uint64_t maskedState;

    //  First commute the k4 operator through the state. We can commute through
    //  until we reach the orbital labelled by k4
    
    //  Mask all occupations to the right of k4
    maskedState = state & (k4Occupied-1);
    
    #if _DEBUG_
        std::cout<<"\n\tMASKED STATE 1: "<<maskedState<<std::endl;
    #endif
    
    //  e.g. if the state is 101001010 and the k4 orbital (in 0,1,2,3,...) is number 
    //  3 then we mask with  000000111 and obtain 000000010.  
	
    //  Then simply count the number of of non-zero bits in that result and see 
    //  if the number is even or odd at the end
    
    #if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient = utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient = utilities::binary::HammingWeight64(maskedState);
    
    #endif

    //  Remove the k4 orbital from the state (with bitwise XOR), 
	//	since it's been annihilated
    outState = state ^ k4Occupied;

	//  Next we commute the k3 operator through the remaining state
    maskedState = outState & (k3Occupied-1);
    
    #if _DEBUG_
        std::cout<<"\n\tMASKED STATE 2: "<<maskedState<<std::endl;
    #endif
    
    #if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient += utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient += utilities::binary::HammingWeight64(maskedState);
    
    #endif

	//  Annihilate the k3 orbital to produce the output state (with bitwise XOR)
	
	outState = outState ^ k3Occupied;	
	
	//	We generate the final state by adding in the k1 and k2 labelled
	//	creation operators. We add k2 first (because it's to the right of k1), 
	//	then k1 afterwards. The procedure is (almost) identical to the one for
	//	k3 and k4 used above

	//  Mask all occupations to the right of k2 
    maskedState = outState & (k2Occupied-1);
    
    #if _DEBUG_
        std::cout<<"\n\tMASKED STATE 3: "<<maskedState<<std::endl;
    #endif

	//	Count the number of bits in the remaining state
	#if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient += utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient += utilities::binary::HammingWeight64(maskedState);
    
    #endif

	//	Add the k2 orbital to the state (with bitwise OR)
	outState = outState | k2Occupied;
	
	//	Finally, add in k1
	
	//  Mask all occupations to the right of k1
    maskedState = outState & (k1Occupied-1);
    
    #if _DEBUG_
        std::cout<<"\n\tMASKED STATE 4: "<<maskedState<<std::endl;
    #endif

	//	Count the number of bits in the remaining state
	#if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient += utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient += utilities::binary::HammingWeight64(maskedState);
    
    #endif
    
    #if _DEBUG_
        std::cout<<"\n\tcoefficient = "<<coefficient<<std::endl;
    #endif
	
	//	Add the k1 orbital to the state (with bitwise OR)
	outState = outState | k1Occupied;

    //  Convert our coefficient into a + or - sign only. 
    if(coefficient & 1)
    {
        coefficient = -1;
    }
    else
    {
        coefficient = +1;
    }

    return  outState;
}
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Apply an operator c^+_k1 c_k2 to a given state 
//!
//!	NOTE: the states are assumed to be written in ascending order of orbital labels
//!
//! NOTE: this implementation assumes a 64-bit architecture
//!
//! \return The index of the final state.
////////////////////////////////////////////////////////////////////////////////

inline uint64_t Evaluate_CdC(
    int&  coefficient,              //!<    Address to store the sign of the output state  
    const uint64_t k1Occupied,      //!<    First orbital occupied 
    const uint64_t k2Occupied,      //!<    Second orbital occupied   
    const uint64_t state,           //!<    The state to be operated on
    const uint64_t highestOrbital)  //!<    Highest possible orbital
{
	//	Store the following results for subsequent use. These values store
	//	e.g. ...00000010000...
	//				  ^ - corresponds to the k1 or k2 orbital 
	//	Note that we denote ...0001 as the 0th orbital

    #if _DEBUG_
        std::cout<<"k1: "<<std::bitset<64>(k1Occupied)<<" k2: "<<std::bitset<64>(k2Occupied)<<" state: "<<state<<std::endl;
    #endif

	//  Check conditions for which the result will be zero
    //  NOTE: lowest state is const static uint64_t 1
    if( 
    (k2Occupied>highestOrbital)  //  Orbitals higher than highestOrbital are not occupied                                               
    || ((state & k2Occupied)==0) //  If the state does not already contain k2
							     //  then c_k2 acts to give zero 
	|| (!(k1Occupied & k2Occupied) && (state & k1Occupied)!=0) //	 If the c^+_k1 acts on an identical
								 //  state then we get zero, unless k1 and k2 are equal
    )
    {   
        coefficient = 0;
        
        #if _DEBUG_
            std::cout<<"return: 0"<<std::endl;
        #endif
        
        return 0;    //  Default return, but it won't be used   
    }
	
	//  If the input state is not zero then the output state is related to 
    //  it by a simple set of bitshift operations.
    //  We also need to keep track of the sign of the output state

    if(k1Occupied & k2Occupied)
	{
		//	In this case the output state is the same as the input state
		//	and the sign is always positive
		
		coefficient = 1;
		
		#if _DEBUG_
		    std::cout<<"return: 1"<<std::endl;
		#endif
		
		return state;
	}

	//	Declare some working variables
	
	uint64_t outState;
	uint64_t maskedState;

    //  First commute the k2 operator through the state. We can commute through
    //  until we reach the orbital labelled by k2
    
    //  Mask all occupations to the right of k2
    maskedState = state & (k2Occupied-1);
	
	//  Then simply count the number of of non-zero bits in that result and see 
    //  if the number is even or odd at the end
	#if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient = utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient = utilities::binary::HammingWeight64(maskedState);
    
    #endif
	
	//  Remove the k2 orbital from the state (with bitwise XOR), 
	//	since it's been annihilated
    outState = state ^ k2Occupied;
	
	//	We generate the final state by adding in the k1 labelled
	//	creation operator. 
	
	//  Mask all occupations to the right of k1 
    maskedState = outState & (k1Occupied-1);

	//	Count the number of bits in the remaining state
	#if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient += utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient += utilities::binary::HammingWeight64(maskedState);
    
    #endif
    
    #if _DEBUG_
        std::cout<<"\n\tcoefficient = "<<coefficient<<std::endl;
    #endif
	
	//	Add the k1 orbital to the state (with bitwise OR)
	outState = outState | k1Occupied;
	
	 //  Convert our coefficient into a + or - sign only. 
    if(coefficient & 1)
    {
        coefficient = -1;
    }
    else
    {
        coefficient = +1;
    }

	return outState;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Apply an operator c_k1 c^+_k2 to a given state 
//!
//!	NOTE: the states are assumed to be written in ascending order of orbital labels
//!
//! NOTE: this implementation assumes a 64-bit architecture
//!
//! \return The index of the final state.
////////////////////////////////////////////////////////////////////////////////

inline uint64_t Evaluate_CCd(
    int&  coefficient,              //!<    Address to store the sign of the output state  
    const uint64_t k1Occupied,      //!<    First orbital occupied 
    const uint64_t k2Occupied,      //!<    Second orbital occupied   
    const uint64_t state,           //!<    The state to be operated on
    const uint64_t highestOrbital)  //!<    Highest possible orbital
{
	//	Store the following results for subsequent use. These values store
	//	e.g. ...00000010000...
	//				  ^ - corresponds to the k1 or k2 orbital 
	//	Note that we denote ...0001 as the 0th orbital

    #if _DEBUG_
        std::cout<<"k1: "<<std::bitset<64>(k1Occupied)<<" k2: "<<std::bitset<64>(k2Occupied)<<" state: "<<state<<std::endl;
    #endif

	//  Check conditions for which the result will be zero
    //  NOTE: lowest state is const static uint64_t 1
    if( 
    (k2Occupied>highestOrbital)  //  Orbitals higher than highestOrbital are not occupied                                               
    || ((state & k2Occupied)!=0) //  If the state already contains k2
							     //  then c^+_k2 acts to give zero 
    )
    {   
        coefficient = 0;
        
        #if _DEBUG_
            std::cout<<"return: 0"<<std::endl;
        #endif
        
        return 0;    //  Default return, but it won't be used   
    }
	
	//  If the input state is not zero then the output state is related to 
    //  it by a simple set of bitshift operations.
    //  We also need to keep track of the sign of the output state

    if(k1Occupied & k2Occupied)
	{
		//	In this case the output state is the same as the input state
		//	and the sign is always positive
		
		coefficient = 1;
		
		#if _DEBUG_
		    std::cout<<"return: 1"<<std::endl;
		#endif
		
		return state;
	}

	//	Declare some working variables
	
	uint64_t outState = state;
	uint64_t maskedState;
    
    //  First add the k2 orbital to the state
    
    //  Mask all occupations to the right of k2 
    maskedState = outState & (k2Occupied-1);

	//	Count the number of bits in the remaining state
	#if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient = utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient = utilities::binary::HammingWeight64(maskedState);
    
    #endif
    
    #if _DEBUG_
        std::cout<<"\n\tcoefficient = "<<coefficient<<std::endl;
    #endif
	
	//	Add the k2 orbital to the state (with bitwise OR)
	outState = outState | k2Occupied;
    
    //  Next commute the k1 operator through the state. We can commute through
    //  until we reach the orbital labelled by k1
    
    //  Mask all occupations to the right of k1
    maskedState = outState & (k1Occupied-1);
	
	//  Then simply count the number of of non-zero bits in that result and see 
    //  if the number is even or odd at the end
	#if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
    
        coefficient += utilities::binary::HammingWeight64Iterative(maskedState);
        
    #else
    
        coefficient += utilities::binary::HammingWeight64(maskedState);
    
    #endif
	
	//  Remove the k1 orbital from the state (with bitwise XOR), 
	//	since it's been annihilated
    outState ^= k1Occupied;

	 //  Convert our coefficient into a + or - sign only. 
    if(coefficient & 1)
    {
        coefficient = -1;
    }
    else
    {
        coefficient = +1;
    }

	return outState;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

};  //  End namespace utilities 		
			
#endif		
