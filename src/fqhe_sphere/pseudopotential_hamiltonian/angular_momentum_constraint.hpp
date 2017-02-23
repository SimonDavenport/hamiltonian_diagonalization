////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport       
//!                                                                             
//!	 \file
//!     This file defines a function to impose a particular angular momentum 
//!     sector for Fock states in the FQHE sphere geometry
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

#ifndef _ANGULAR_MOMENTUM_CONSTRAINTS_HPP_INCLUDED_
#define _ANGULAR_MOMENTUM_CONSTRAINTS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/mathematics/binary_number_tools.hpp" 

namespace diagonalization
{ 
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to test whether a given Fock state is in a given  
    //! quantum number sector.
    //!
    //! This function defines the representation of the Fock states in terms of
    //! angular momentum quantum numbers m
    ////////////////////////////////////////////////////////////////////////////////
    inline bool TestAngularMomentumSector(
        const fock_t state,         //!<    Fock state to be tested
        const mState_t lzTot,       //!<    Given angular momentum sector
                                    //!     to test against
        const mState_t max)         //!<    Bounding angular momentum quantum number          
    {
        //  This function assumes that the Fock states are labelled by
        //  angular momentum quantum numbers from -m to m in steps of 2
        fock_t tempState = state;
        //  Iterate over occupied states to find the total angular momentum value
        mState_t cumulativeLz = 0;
        while(0 != tempState)
        {
            //  Isolate right-most bit and determine its position
            fock_t rightMost = tempState & -tempState;
            cumulativeLz += 2*utilities::binary::Base2Log(rightMost) - max;  
            //  Mask off the current right-most bit and iterate
            tempState = tempState ^ rightMost;
        }
        //  Compare with the actual lzTot sector
        return (lzTot == cumulativeLz);
    }
}   //  End namespace diagonalization
#endif
