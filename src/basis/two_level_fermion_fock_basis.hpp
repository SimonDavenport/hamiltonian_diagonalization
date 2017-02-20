////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 03/02/2015
//!
//!  \file
//!		This file defines a container for a Fock basis of fermions in a 
//!     two level system
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

#ifndef _TWO_LEVEL_FERMION_FOCK_BASIS_HPP_INCLUDED_
#define _TWO_LEVEL_FERMION_FOCK_BASIS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../utilities/general/orbital_and_state_defs.hpp"
                                                //  Define fock_t and fock_t
#include "../utilities/wrappers/mpi_wrapper.hpp"        
                                                //  For mpi functionality
#include "../utilities/algorithms/binary_search.hpp"      
                                                //  For binary search algorithm
#include "../utilities/mathematics/binomial_table.hpp"     
                                                //  For binomials
#include "../utilities/mathematics/binary_number_tools.hpp"
                                                //  For bit counting and bit 
                                                //  permutation operations
#include <boost/function.hpp>                   //  For passing quantum number constraint function
#include <boost/bind.hpp>                       //  For quantum number constraint parameters
#include <string.h>                             //  For memcpy function
#include <algorithm>                            //  For std::sort

#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{

////////////////////////////////////////////////////////////////////////////////
//! \brief A struct to store a pair of binary numbers to indicate the
//! occupations in a 2 level system
//!
////////////////////////////////////////////////////////////////////////////////

struct TwoLevelFockState
{
    fock_t m_first;
    fock_t m_second;

    //  Default constructor
    TwoLevelFockState();

    //  Constructor from two separate values  
    TwoLevelFockState(const fock_t first,const fock_t second);

    //  Comparison operator (required for binary search)
    friend bool operator > (const TwoLevelFockState& state1,const TwoLevelFockState& state2);
    
    //  Comparison operator (required for binary search)
    friend bool operator < (const TwoLevelFockState& state1,const TwoLevelFockState& state2);
    
    //  Equals operator (required for binary search)
    friend bool operator == (const TwoLevelFockState& state1,const TwoLevelFockState& state2);
    
    //  Overload out stream operator
    friend std::ostream& operator<< (std::ostream& os, const TwoLevelFockState& state);
};

////////////////////////////////////////////////////////////////////////////////
//! \brief Store the Fermion Fock states as two 64 bit integer arrays
//!
////////////////////////////////////////////////////////////////////////////////

class TwoLevelFermionFockBasis
{
    public:

    std::vector<TwoLevelFockState> m_basis;  
                                        //!<    Array of two level Fock states
    bool m_basisStored;                 //!<    Set to true if a Fock basis is
                                        //!     stored
                                        
    fock_t m_firstState;                //!<    Index of the first state stored
                                        //!     (used for parallelized implementations)                               
    public:
    
    //  Constructors
    TwoLevelFermionFockBasis();
    TwoLevelFermionFockBasis(const TwoLevelFermionFockBasis& other);
    TwoLevelFermionFockBasis& operator=(const TwoLevelFermionFockBasis& other);
    
    //  Destructor
    ~TwoLevelFermionFockBasis();

    void Clear();
    
    fock_t GetFullFockSpaceDimension(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals1,const iSize_t nbrOrbitals2) const;
    
    fock_t GetReducedFockSpaceDimension(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals1,const iSize_t nbrOrbitals2,const utilities::MpiWrapper& mpi) const;

    void AddState(const fock_t state1,const fock_t state2);
   
    void FindFockStateIndex(const fock_t offset,const TwoLevelFockState& state,
        fock_t& ) const;
    
    void GetFirstState(const fock_t firstNodeIndex,
        std::vector<TwoLevelFockState>::const_iterator& it_state) const;
    
    void GetNextState(std::vector<TwoLevelFockState>::const_iterator& it_state) const;
    
    bool GenerateFockSpace(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals1,const iSize_t nbrOrbitals2,
        const std::function<bool(const TwoLevelFockState& state)>& Constraint,utilities::MpiWrapper& mpi);
    
    fock_t GetStoredDimension() const;
    
    void GetFockBasis(fock_t* buffer1,fock_t* buffer2,const fock_t dim,utilities::MpiWrapper& mpi) const;
    
    void SetFockBasis(fock_t* buffer1,fock_t* buffer2,const fock_t dim,const fock_t firstState);
    
};

}   //  End diagonalization namespace

#endif
