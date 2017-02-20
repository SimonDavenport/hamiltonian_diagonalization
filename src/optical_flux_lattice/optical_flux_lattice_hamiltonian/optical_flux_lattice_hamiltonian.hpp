////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 02/02/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store the optical flux lattice interacting
//!     Hamiltonian. See e.g. PRL 109, 265301(2012) 	   
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

#ifndef _OPTICAL_FLUX_LATTICE_HAMILTONIAN_HPP_INCLUDED_
#define _OPTICAL_FLUX_LATTICE_HAMILTONIAN_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <fstream>                      //  for fstream
#include <vector>                       //  for std::vector
#include <boost/program_options.hpp>    //  for program option handling
#include <functional>                   //  for passing quantum number constraint function
#include <string.h>                     //  for memcpy function

//  Include declaration of associated program options
#include "../program_options/optical_flux_lattice_hamiltonian_options.hpp"
#include "../program_options/single_particle_hamiltonian_options.hpp"

//  Single particle model definition
#include "../single_particle_hamiltonian/single_particle_hamiltonian.hpp"
#include "../single_particle_hamiltonian/single_particle_hamiltonian_array.hpp"

//  Interacting Hamiltonian parameters
#include "optical_flux_lattice_hamiltonian_data.hpp"

//  Functions to calculate the Hamiltonian
#include "../../hamiltonians/spinless_fermion_hamiltonian.hpp"

//  Functions to calculate observables
#include "../../observables/fermion_observables.hpp"

//  Data structures to store Hamiltonian coefficient data
#include "lookup_tables.hpp"

//  2D Linear momentum constraints functions
#include "linear_momentum_constraints_2D.hpp"

//  Declaration of kState_t, iSize_t and fock_t
#include "../../utilities/general/orbital_and_state_defs.hpp"

//  Matrix-vector routines
#include "../../hamiltonians/matrix_vector_routines.hpp"

//  Include various utility classes and definitions
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/mathematics/lattice_vector.hpp"
#include "../../utilities/general/dcmplx_type_def.hpp"
#include "../../utilities/general/pi_const_def.hpp"
#include "../../utilities/mathematics/array_convolution.hpp"
#include "../../utilities/general/load_bar.hpp"
#include "../../utilities/wrappers/sqlite_wrapper.hpp"
#include "../../utilities/mathematics/binary_number_tools.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../../utilities/mathematics/kahan_arithmetic.hpp"
#include "../../utilities/general/template_tools.hpp"
#include "../../utilities/mathematics/modulo.hpp"

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#include <bitset>
#endif

namespace diagonalization
{

//////      STRUCT/CLASS DECLARATIONS       ////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief The OpticalFluxLatticeHamiltonian class contains a matrix representation of
//! the interacting optical flux lattice Hamiltonian and functions to construct
//! it in a given basis and to diagonalize it. 
//!
//////////////////////////////////////////////////////////////////////////////////

class OpticalFluxLatticeHamiltonian
{
    //  Friend classes
    
    friend class ObservablesManager;

    private:

    //////      Private data structures     ////////////////////////////////////

    SpinlessFermionHamiltonian<dcmplx> m_hamiltonian;   
                                        //!<  An object to contain the Hamiltonian itself
    
    OpticalFluxLatticeHamiltonianData m_params;
                                        //!<  Class internal implementation of the     
                                        //!   OpticalFluxLatticeHamiltonianData struct          
    SingleParticleParameters m_singleParticleData;
                                        //!<  Instance of the single particle
                                        //!   Hamiltonian data struct

    QuadraticLookUpTables m_quadraticTables;    
                                        //!<  Class container for regular array look-up tables
                                        //!   used in the case where we have full 2D momentum conservation 
    QuarticLookUpTables m_quarticTables;        
                                        //!<  Class container for regular array look-up tables
                                        //!   used in the case where we have full 2D momentum conservation 
    
    QuadraticLookUpHashTables m_quadraticHashTables;
                                        //!<  Class container for multi-hash maps:
                                        //!   used when we have only ky-momentum conservation
    QuarticLookUpHashTables m_quarticHashTables;
                                        //!<  Class container for multi-hash maps:
                                        //!   used when we have only ky-momentum conservation

    std::vector<dcmplx> m_basisChangeMatrix;        
                                        //!<    Change of basis matrix;
    
    std::vector<double> m_magnetization;//!<  A list of the magnetization values
    
    //////      Private functions       ////////////////////////////////////////

    //  Generate a table containing each k1 value for a given k2
    void Generate2KTable(utilities::MultiHashMultiMap<kState_t>* kHashTable,const utilities::MpiWrapper& mpi);

    //  Generate an array containing each k1 value for a given k2,k3,k4
    void Generate4KTable(std::vector<kState_t>* kTable,std::vector<int>* gxTable,std::vector<int>* gyTable,const utilities::MpiWrapper& mpi);    

    void Generate4KTable(std::vector<kState_t>* kTable,utilities::MultiHashMultiMap<kState_t>* kHashTable,const utilities::MpiWrapper& mpi);

    //  Find a k1 value for a given k2,k3,k4
    kState_t FindK1(const kState_t k2,const kState_t k3,const kState_t k4,int& gTot1,int& gTot2) const;

    //  Find a k1y value for a given k2y,k3y,k4y assuming only ky conservation
    kState_t FindK1(const kState_t k2y,const kState_t k3y,const kState_t k4y,int& gTot2) const;

    //  Generate Vkkkk from the single particle Bloch coefficients
    void GenerateQuarticTerms(std::vector<dcmplx>* table,const iSize_t dimension,SingleParticleHamiltonian* blochTable,std::vector<int>* gxTable,std::vector<int>* gyTable,utilities::MpiWrapper& mpi) const;

    //  Convert regular array type tables to hash tables
    void ConvertTableFormat(const utilities::MpiWrapper& mpi);
    
	//	Read/write Vkkkk table data to/from files to save re-calculation time
	void TablesToFile(utilities::MpiWrapper& mpi);
	void TablesFromFile(utilities::MpiWrapper& mpi);

	//  Change to a Wannier wave function basis
    void ChangeToWannierBasis(SingleParticleHamiltonian* blochTable,const double matrixElementTol,utilities::MpiWrapper& mpi);

	//  Calculate the sigma^z map and store in m_magnetization
	void StoreSigmaZMap(SingleParticleHamiltonian* blochTable,const utilities::MpiWrapper& mpi);

    //  Calculate and store a list of occupation probabilities
    bool CalculateOccupations(utilities::MpiWrapper& mpi) const;
    
    //  Calculate and store the susceptibility
    bool CalculateSusceptibility(utilities::MpiWrapper& mpi) const;

    //  Generate a plot of the Hamiltonian
    bool PlotHamiltonian(const utilities::MpiWrapper& mpi) const;

    //  Output a list of the most probable Fock states 
    bool ListMostProbable(const iSize_t nbrStates,utilities::MpiWrapper& mpi) const;

    //  Calculate and store the density-density function
    bool CalculateDensityDensityFunction(utilities::MpiWrapper& mpi) const;

    //  Calculate the inverse participation ratio
    bool CalculateParticipationRatio(utilities::MpiWrapper& mpi) const;

    //  Calculate a translational density-density function
    bool CalculateTranslationalDensityDensityFunction(utilities::MpiWrapper& mpi) const;
    
    //  Calculate a rotational density-density function
    bool CalculateRotationalDensityDensityFunction(utilities::MpiWrapper& mpi) const;

    //////      Public interface functions     ////////////////////////////////

    public:
    
    //  Default constructor 
    OpticalFluxLatticeHamiltonian();
    
    //  Constructor   
    OpticalFluxLatticeHamiltonian(OpticalFluxLatticeHamiltonianData params);
    
    //  Constructor from command line arguments
    OpticalFluxLatticeHamiltonian(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);

    //  Destructor
    ~OpticalFluxLatticeHamiltonian();
    
    //  Update SQL status
    void UpdateSqlStatus(const utilities::MpiWrapper& mpi);
    
    //  Build look-up tables
    void BuildLookUpTables(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);
    
    //  Set the sector to be diagonalized
    void SetSector(const kState_t kxTot,const kState_t kyTot);
    
    //  Build Fock basis
    void BuildFockBasis(utilities::MpiWrapper& mpi);
    
    //  Put entries into full Hamiltonian
    void BuildHamiltonian(utilities::MpiWrapper& mpi);

    //  Call diagonalization routine
    void Diagonalize(utilities::MpiWrapper& mpi);

    //  Store the full Hamiltonian in a file (including sparse format)
    void StoreHamiltonian(utilities::MpiWrapper& mpi);
    
    //  Retrieve the full Hamiltonian from a file (including sparse format)
    void RetrieveHamiltonian(utilities::MpiWrapper& mpi);

    //  Destroy the currently stored Hamiltonian
    void ClearHamiltonian();
    
    //  Save/retrieve the current eigensystem to/form a file 
    void EigensystemToFile(const bool writeEigenvalues,const bool writeEigenvectors,utilities::MpiWrapper& mpi);    
    void EigensystemFromFile(const bool readEigenvalues,const bool readEigenvectors,utilities::MpiWrapper& mpi);
};

}   //  End diagonalization namespace

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#endif
