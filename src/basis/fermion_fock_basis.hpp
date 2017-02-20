////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 26/01/2015
//!
//!  \file
//!		This file defines a container for a Fock basis of fermions
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

#ifndef _FERMION_FOCK_BASIS_HPP_INCLUDED_
#define _FERMION_FOCK_BASIS_HPP_INCLUDED_

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
#include <functional>                       //  For passing quantum number constraint function
#include <string.h>                         //  For memcpy function
#include <algorithm>                        //  For std::sort

#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{

////////////////////////////////////////////////////////////////////////////////
//! \brief Store the Fermion Fock states as a 64 bit integer array
//!
////////////////////////////////////////////////////////////////////////////////

class FermionFockBasis
{
    public:

    std::vector<fock_t> m_basis;    //!<    Array of Fock states
    
    bool m_basisStored;             //!<    Set to true if a Fock basis is
                                    //!     stored
                                        
    fock_t m_firstState;            //!<    Index of the first state stored
                                    //!     (used for parallelized implementations)
  
    public:
    
    //  Constructors
    FermionFockBasis();
    FermionFockBasis(const FermionFockBasis& other);
    FermionFockBasis& operator=(const FermionFockBasis& other);
    
    //  Destructor
    ~FermionFockBasis();

    void Clear();
    
    fock_t GetFullFockSpaceDimension(const iSize_t nbrOrbitals,
        const iSize_t nbrParticles) const;
    
    fock_t GetReducedFockSpaceDimension(const iSize_t nbrOrbitals,
        const iSize_t nbrParticles,const utilities::MpiWrapper& mpi) const;

    fock_t FindFockStateIndex(const fock_t offset,
        const fock_t outState) const;
    
    fock_t GetFirstState(const fock_t firstNodeIndex,
        const iSize_t nbrParticles,std::vector<fock_t>::const_iterator& it_state) const;
    
    fock_t GetNextState(const fock_t state,std::vector<fock_t>::const_iterator& it_state) const;
    
    bool GenerateFockSpace(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals,
        const std::function<bool(const fock_t state)>& Constraint,utilities::MpiWrapper& mpi);
    
    fock_t GetStoredDimension() const;
    
    void GetFockBasis(fock_t* buffer,const fock_t dim,utilities::MpiWrapper& mpi) const;
    
    void SetFockBasis(fock_t* buffer,const fock_t dim,const fock_t firstState);
        
};

}   //  End namespace diagonalization

#endif
