////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                           
//!                                                                             
//!                      \date Last Modified: 08/06/2015 
//!                                                                             
//!	 \file
//!     This file defines functions that impose 2D total linear momentum
//!     constraints on a given fock state
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "linear_momentum_constraints_2D.hpp"

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief A function to test whether a given Fock state is in a given quantum 
//! number sector.
//!
//! This function defines the representation of the Fock states in terms of
//! linear momentum quantum numbers kx,ky.
//!
////////////////////////////////////////////////////////////////////////////////

bool TestLinearMomentumSector2(
    const fock_t state,         //!<    Fock state to be tested
    const kState_t kxTot,       //!<    Given linear momentum sector
                                //!<    To test against
    const kState_t kyTot,       //!<    Given linear momentum sector
                                //!<    To test against
    const kState_t dimX,        //!<    X dimension of k-space grid
    const kState_t dimY)        //!<    Y dimension of k-space grid                                    
{
    //  The occupation basis is of 64-bit integer format e.g. 000100101...
    //  
    //  Let's say that the ith bit from the right corresponds to the 
    //  kx*ny+ky th state, where ny is the number of ky states.
    //  e.g. the i=3 state (with ny=4 say) corresponds to kx=0,ky=3
    //  the i=5 state would in that case be kx=1,ky=0 etc.
    //
    //  States in a given sector have the same total value of kx and
    //  ky (modulo nx or ny respectively)
    //
    //  The number of sectors is equal to nx*ny. We shall label different
    //  sectors by a single quantum number kxTot*ny + kyTot
    
    //  Determine the kx_tot and ky_tot for the specified state

    fock_t tempState = state;

    int x = 0;
    int y = 0;

    while(0 != tempState)
    {
        //  Isloate right-most bit and determine its position
        
        fock_t rightMost = tempState & -tempState;

        kState_t position = utilities::binary::Base2Log(rightMost);

        //  Convert to linear momentum quantum numbers
        //  and accumulate total linear momentum

        x += floor((double)position/(dimY));
        y += utilities::Modulo<int>(position,dimY);
 
        //  Mask off the current right-most bit and iterate
        
        tempState = tempState ^ rightMost;
    }
    
    x = utilities::Modulo<int>(x,dimX);
    y = utilities::Modulo<int>(y,dimY);
    
    x = x<0 ? x + dimX : x;
    y = y<0 ? y + dimY : y;
    
    //  Compare the specified sector with the given linear momentum
    //  sector

    return (x == kxTot && y == kyTot);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief A function to test whether a given Fock state is in a given quantum 
//! number sector.
//!
//! This function defines the representation of the Fock states in terms of
//! linear momentum quantum number ky only.
//!
////////////////////////////////////////////////////////////////////////////////

bool TestLinearMomentumSector1(
    const fock_t state,         //!<    Fock state to be tested
    const kState_t kyTot,       //!<    Given linear momentum sector
                                //!<    To test against
    const kState_t dimY)        //!<    Y dimension of k-space grid                                    
{
    //  The occupation basis is of 64-bit integer format e.g. 000100101...
    //  
    //  Let's say that the ith bit from the right corresponds to the 
    //  kx*ny+ky th state, where ny is the number of ky states.
    //  e.g. the i=3 state (with ny=4 say) corresponds to kx=0,ky=3
    //  the i=5 state would in that case be kx=1,ky=0 etc.
    //
    //  States in a given sector have the same total value of kx and
    //  ky (modulo nx or ny respectively)
    //
    //  The number of sectors is equal to nx*ny. We shall label different
    //  sectors by a single quantum number kxTot*ny + kyTot
    
    //  Determine the ky_tot for the specified state

    fock_t tempState = state;

    int y = 0;

    while(0 != tempState)
    {
        //  Isloate right-most bit and determine its position
        
        fock_t rightMost = tempState & -tempState;

        kState_t position = utilities::binary::Base2Log(rightMost);

        //  Convert to linear momentum quantum numbers
        //  and accumulate total linear momentum

        y += utilities::Modulo<int>(position,dimY);
 
        //  Mask off the current right-most bit and iterate
        
        tempState = tempState ^ rightMost;
    }
    
    y = utilities::Modulo<int>(y,dimY);
    
    y = y<0 ? y + dimY : y;
    
    //  Compare the specified sector with the given linear momentum
    //  sector

    return y == kyTot;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
////////////////////////////////////////////////////////////////////////////////
//! \brief A function to add two k-states to give a result that is periodic
//! in the unit cell
//!
//! \return k1+k2 pariodic in 2D momentum space
//!
////////////////////////////////////////////////////////////////////////////////

kState_t StateAddition(
    const kState_t k1,          //!<    First k state
    const kState_t k2,          //!<    First k state
    const kState_t dimX1,       //!<    X dimension of k-space grid
    const kState_t dimY1,       //!<    Y dimension of k-space grid
    const kState_t dimX2,       //!<    X dimension of k-space grid
    const kState_t dimY2)       //!<    Y dimension of k-space grid
{
    //  Extract the kx and ky components of k1 and k2
    
    int k1x = floor((double)k1/(dimY1));
    int k1y = utilities::Modulo<int>(k1,dimY1);
    
    int k2x = floor((double)k2/(dimY2));
    int k2y = utilities::Modulo<int>(k2,dimY2);
    
    //  Perform modulo addition of the x and y components
    
    int k3x = utilities::Modulo<int>(k1x+k2x,dimX1);
    int k3y = utilities::Modulo<int>(k1y+k2y,dimY1);
    
    //  Recombine to an output state
    
    return (kState_t)(k3x*dimY1+k3y);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
////////////////////////////////////////////////////////////////////////////////
//! \brief A function to subtract two k-states to give a result that is periodic
//! in the unit cell
//!
//! \return k1-k2 pariodic in 2D momentum space
//!
////////////////////////////////////////////////////////////////////////////////

kState_t StateSubtraction(
    const kState_t k1,          //!<    First k state
    const kState_t k2,          //!<    First k state
    const kState_t dimX1,       //!<    X dimension of k-space grid
    const kState_t dimY1,       //!<    Y dimension of k-space grid
    const kState_t dimX2,       //!<    X dimension of k-space grid
    const kState_t dimY2)       //!<    Y dimension of k-space grid
{
    //  Extract the kx and ky components of k1 and k2
    
    int k1x = floor((double)k1/(dimY1));
    int k1y = utilities::Modulo<int>(k1,dimY1);
    
    int k2x = floor((double)k2/(dimY2));
    int k2y = utilities::Modulo<int>(k2,dimY2);
    
    //  Perform modulo addition of the x and y components
    
    int k3x = utilities::Modulo<int>(k1x-k2x,dimX1);
    int k3y = utilities::Modulo<int>(k1y-k2y,dimY1);
    
    k3x = k3x<0 ? k3x + dimX1 : k3x;
    k3y = k3y<0 ? k3y + dimY1 : k3y;
    
    //  Recombine to an output state
    
    return (kState_t)(k3x*dimY1+k3y);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
////////////////////////////////////////////////////////////////////////////////
//! \brief A function to impose a number of R60 rotations on a given k state
//!
//! \return R60 operator applied r times to the state
//!
////////////////////////////////////////////////////////////////////////////////

kState_t R60(
    const kState_t k,           //!<    Input k state
    const int r,                //!<    Number of times to apply the rotation
    const kState_t dimX,        //!<    X dimension of k-space grid
    const kState_t dimY)        //!<    Y dimension of k-space grid
{
    //  Extract the kx and ky components of k
    
    int kx = floor((double)k/(dimY));
    int ky = utilities::Modulo<int>(k,dimY);

    int k1x = 0;
    int k1y = 0;

    switch(utilities::Modulo<int>(r,6))
    {
       case 0:
            k1x = kx;
            k1y = ky;
            break;
       case 1:
            k1x = utilities::Modulo<int>(kx+ky,dimX);
            k1y = utilities::Modulo<int>(-kx,dimY);
            break;
       case 2:
            k1x = utilities::Modulo<int>(ky,dimX);
            k1y = utilities::Modulo<int>(-kx-ky,dimY);
            break;
       case 3:
            k1x = utilities::Modulo<int>(-kx,dimX);
            k1y = utilities::Modulo<int>(-ky,dimY);
            break;
       case 4:
            k1x = utilities::Modulo<int>(-kx-ky,dimX);
            k1y = utilities::Modulo<int>(kx,dimY);
            break;
       case 5:
            k1x = utilities::Modulo<int>(-ky,dimX);
            k1y = utilities::Modulo<int>(kx+ky,dimY);
            break;
    }
    
    return (kState_t)(k1x*dimY+k1y);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
////////////////////////////////////////////////////////////////////////////////
//! \brief A function to impose a number of inverse R60 rotations on a given 
//! k state
//!
//! \return Inverse R60 operator applied r times to the state
//!
////////////////////////////////////////////////////////////////////////////////

kState_t InverseR60(
    const kState_t k,           //!<    Input k state
    const int r,                //!<    Number of times to apply the rotation
    const kState_t dimX,        //!<    X dimension of k-space grid
    const kState_t dimY)        //!<    Y dimension of k-space grid
{
    //  Extract the kx and ky components of k
    
    int kx = floor((double)k/(dimY));
    int ky = utilities::Modulo<int>(k,dimY);

    int k1x = 0;
    int k1y = 0;

    switch(utilities::Modulo<int>(r,6))
    {
       case 0:
            k1x = kx;
            k1y = ky;
            break;
       case 1:
            k1x = utilities::Modulo<int>(-ky,dimX);
            k1y = utilities::Modulo<int>(kx+ky,dimY);
            break;
       case 2:
            k1x = utilities::Modulo<int>(-kx-ky,dimX);
            k1y = utilities::Modulo<int>(kx,dimY);
            break;
       case 3:
            k1x = utilities::Modulo<int>(-kx,dimX);
            k1y = utilities::Modulo<int>(-ky,dimY);
            break;
       case 4:
            k1x = utilities::Modulo<int>(ky,dimX);
            k1y = utilities::Modulo<int>(-kx-ky,dimY);
            break;
       case 5:
            k1x = utilities::Modulo<int>(kx+ky,dimX);
            k1y = utilities::Modulo<int>(-kx,dimY);
            break;
    }
    
    return (kState_t)(k1x*dimY+k1y);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace diagonalization
