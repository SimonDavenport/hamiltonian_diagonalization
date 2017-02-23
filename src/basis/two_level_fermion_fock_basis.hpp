////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file defines a container for a Fock basis of fermions in a 
//!     two level system
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

#ifndef _TWO_LEVEL_FERMION_FOCK_BASIS_HPP_INCLUDED_
#define _TWO_LEVEL_FERMION_FOCK_BASIS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../utilities/general/orbital_and_state_defs.hpp"
#include "../utilities/wrappers/mpi_wrapper.hpp"
#include "../utilities/algorithms/binary_search.hpp"
#include "../utilities/mathematics/binomial_table.hpp"
#include "../utilities/mathematics/binary_number_tools.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp> 
#include <string.h> 
#include <algorithm>
#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    //!
    //! A struct to store a pair of binary numbers to indicate the
    //! occupations in a 2 level system
    //!
    struct TwoLevelFockState
    {
        fock_t m_first;
        fock_t m_second;
        TwoLevelFockState();
        TwoLevelFockState(const fock_t first, const fock_t second);
        friend bool operator > (const TwoLevelFockState& state1, const TwoLevelFockState& state2);
        friend bool operator < (const TwoLevelFockState& state1, const TwoLevelFockState& state2);
        friend bool operator == (const TwoLevelFockState& state1, const TwoLevelFockState& state2);
        friend std::ostream& operator<< (std::ostream& os, const TwoLevelFockState& state);
    };

    //!
    //! Store the Fermion Fock states as two 64 bit integer arrays
    //!
    class TwoLevelFermionFockBasis
    {
        public:
        std::vector<TwoLevelFockState> m_basis;  
                                            //!<    Array of two level Fock states
        bool m_basisStored;                 //!<    Set to true if a Fock basis is
                                            //!     stored
                                            
        fock_t m_firstState;                //!<    Index of the first state stored
                                            //!     (used for parallelized implementations)                               
        TwoLevelFermionFockBasis();
        TwoLevelFermionFockBasis(const TwoLevelFermionFockBasis& other);
        TwoLevelFermionFockBasis& operator=(const TwoLevelFermionFockBasis& other);
        ~TwoLevelFermionFockBasis();
        void Clear();
        fock_t GetFullFockSpaceDimension(const iSize_t nbrParticles,
                                         const iSize_t nbrOrbitals1,
                                         const iSize_t nbrOrbitals2) const;
        fock_t GetReducedFockSpaceDimension(const iSize_t nbrParticles,
                                            const iSize_t nbrOrbitals1,
                                            const iSize_t nbrOrbitals2,
                                            const utilities::MpiWrapper& mpi) const;
        void AddState(const fock_t state1, const fock_t state2);
        void FindFockStateIndex(const fock_t offset, const TwoLevelFockState& state,
                                fock_t& ) const;
        void GetFirstState(const fock_t firstNodeIndex,
                           std::vector<TwoLevelFockState>::const_iterator& it_state) const;
        void GetNextState(std::vector<TwoLevelFockState>::const_iterator& it_state) const;
        bool GenerateFockSpace(const iSize_t nbrParticles, const iSize_t nbrOrbitals1,
                               const iSize_t nbrOrbitals2, 
                               const std::function<bool(const TwoLevelFockState& state)>& Constraint, 
                               utilities::MpiWrapper& mpi);
        fock_t GetStoredDimension() const;
        void GetFockBasis(fock_t* buffer1, fock_t* buffer2, const fock_t dim, utilities::MpiWrapper& mpi) const;
        void SetFockBasis(fock_t* buffer1, fock_t* buffer2, const fock_t dim, const fock_t firstState);
    };
}   //  End diagonalization namespace
#endif
