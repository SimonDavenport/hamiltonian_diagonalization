////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 13/02/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store a FQHE Haldane pseduopotential
//!     Hamiltonian for in the sphere geometry. Both 1 and 2LL systems
//!     can be treated.
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

#ifndef _TWO_LEVEL_PSEUDOPOTENTIAL_HAMILTONIAN_HPP_INCLUDED_
#define _TWO_LEVEL_PSEUDOPOTENTIAL_HAMILTONIAN_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

//  Include base class declaration
#include "pseudopotential_hamiltonian_base.hpp"

//  Functions to perform exact diagonalization
#include "../../hamiltonians/two_level_spinless_fermion_hamiltonian.hpp"

//  Data structures to store Hamiltonian coefficient data
#include "lookup_tables.hpp"

//  Angular momentum sector test function
#include "two_level_angular_momentum_constraint.hpp"

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#include <bitset>
#endif

namespace diagonalization
{

//////////////////////////////////////////////////////////////////////////////////
//! \brief The SphereTwoLevelPseudopotentialHamiltonian class defines a FQHE
//! Haldane pseudopotential Hamiltonian the sphere geometry. Both 1 and 2LL
//! systems can be treated.
//!
//////////////////////////////////////////////////////////////////////////////////

class SphereTwoLevelPseudopotentialHamiltonian : 
public SpherePseudopotentialHamiltonianBase<TwoLevelSpinlessFermionHamiltonian<double> >
{
    private:

    //////      Private data structures     ////////////////////////////////////

    QuadraticLookUpTables m_quadraticTables;
                                        //!<  Class container for regular array look-up tables
    QuarticLookUpTables m_quarticTables;//!<  Class container for regular array look-up tables

    QuadraticLookUpTables m_quadraticTables2LL;
                                        //!<  Class container for regular array look-up tables
                                        //!   for 2nd LL interactions
    QuarticLookUpTables m_quarticTables2LL;
                                        //!<  Class container for regular array look-up tables
                                        //!   for 2nd LL interactions
    
    //////      Public interface functions     ////////////////////////////////

    public:

    //  Constructor from minimal data
    SphereTwoLevelPseudopotentialHamiltonian(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals,const std::vector<double>& pseudopotentials,
        const std::vector<double>& pseudopotentials2LL);
    
    //  Constructor from command line arguments
    SphereTwoLevelPseudopotentialHamiltonian(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);

    //  Build look-up tables
    void BuildLookUpTables(const utilities::MpiWrapper& mpi);

    //  Build Fock basis
    void BuildFockBasis(utilities::MpiWrapper& mpi);

    //  Set optional CdC terms
    void SetOccupationEnergies(double* energyLevels,const iSize_t dim,const utilities::MpiWrapper& mpi);
    
    //  Set Fock basis from external array
    void SetFockBasis(fock_t* buffer1,fock_t* buffer2,const fock_t dim);
    
    //  Put entries into full Hamiltonian
    void BuildHamiltonian(utilities::MpiWrapper& mpi);
    
    //  Diagonalize the Hamiltonian
    void Diagonalize(utilities::MpiWrapper& mpi);
};

}   //  End diagonalization namespace

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#endif
