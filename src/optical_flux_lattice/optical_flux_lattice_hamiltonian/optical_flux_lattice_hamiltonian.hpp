////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store the optical flux lattice interacting
//!     Hamiltonian. See e.g. PRL 109, 265301(2012) 	   
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

#ifndef _OPTICAL_FLUX_LATTICE_HAMILTONIAN_HPP_INCLUDED_
#define _OPTICAL_FLUX_LATTICE_HAMILTONIAN_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include <functional>
#include <string.h>
#include "../program_options/optical_flux_lattice_hamiltonian_options.hpp"
#include "../program_options/single_particle_hamiltonian_options.hpp"
#include "../single_particle_hamiltonian/single_particle_hamiltonian.hpp"
#include "../single_particle_hamiltonian/single_particle_hamiltonian_array.hpp"
#include "optical_flux_lattice_hamiltonian_data.hpp"
#include "../../hamiltonians/spinless_fermion_hamiltonian.hpp"
#include "../../observables/fermion_observables.hpp"
#include "lookup_tables.hpp"
#include "linear_momentum_constraints_2D.hpp"
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../hamiltonians/matrix_vector_routines.hpp"
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
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief The OpticalFluxLatticeHamiltonian class contains a matrix representation of
    //! the interacting optical flux lattice Hamiltonian and functions to construct
    //! it in a given basis and to diagonalize it. 
    //////////////////////////////////////////////////////////////////////////////////
    class OpticalFluxLatticeHamiltonian
    {
        friend class ObservablesManager;
        private:
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
        void Generate2KTable(utilities::MultiHashMultiMap<kState_t>* kHashTable, const utilities::MpiWrapper& mpi);
        void Generate4KTable(std::vector<kState_t>* kTable,std::vector<int>* gxTable,std::vector<int>* gyTable, const utilities::MpiWrapper& mpi);    
        void Generate4KTable(std::vector<kState_t>* kTable,utilities::MultiHashMultiMap<kState_t>* kHashTable, const utilities::MpiWrapper& mpi);
        kState_t FindK1(const kState_t k2, const kState_t k3, const kState_t k4, int& gTot1, int& gTot2) const;
        kState_t FindK1(const kState_t k2y, const kState_t k3y, const kState_t k4y, int& gTot2) const;
        void GenerateQuarticTerms(std::vector<dcmplx>* table, const iSize_t dimension,SingleParticleHamiltonian* blochTable,std::vector<int>* gxTable, std::vector<int>* gyTable, utilities::MpiWrapper& mpi) const;
        void ConvertTableFormat(const utilities::MpiWrapper& mpi);
	    void TablesToFile(utilities::MpiWrapper& mpi);
	    void TablesFromFile(utilities::MpiWrapper& mpi);
        void ChangeToWannierBasis(SingleParticleHamiltonian* blochTable, const double matrixElementTol, utilities::MpiWrapper& mpi);
	    void StoreSigmaZMap(SingleParticleHamiltonian* blochTable, const utilities::MpiWrapper& mpi);
        bool CalculateOccupations(utilities::MpiWrapper& mpi) const;
        bool CalculateSusceptibility(utilities::MpiWrapper& mpi) const;
        bool PlotHamiltonian(const utilities::MpiWrapper& mpi) const;
        bool ListMostProbable(const iSize_t nbrStates,utilities::MpiWrapper& mpi) const;
        bool CalculateDensityDensityFunction(utilities::MpiWrapper& mpi) const;
        bool CalculateParticipationRatio(utilities::MpiWrapper& mpi) const;
        bool CalculateTranslationalDensityDensityFunction(utilities::MpiWrapper& mpi) const;
        bool CalculateRotationalDensityDensityFunction(utilities::MpiWrapper& mpi) const;
        public:
        OpticalFluxLatticeHamiltonian();  
        OpticalFluxLatticeHamiltonian(OpticalFluxLatticeHamiltonianData params);
        OpticalFluxLatticeHamiltonian(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        ~OpticalFluxLatticeHamiltonian();
        void UpdateSqlStatus(const utilities::MpiWrapper& mpi);
        void BuildLookUpTables(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void SetSector(const kState_t kxTot,const kState_t kyTot);
        void BuildFockBasis(utilities::MpiWrapper& mpi);
        void BuildHamiltonian(utilities::MpiWrapper& mpi);
        void Diagonalize(utilities::MpiWrapper& mpi);
        void StoreHamiltonian(utilities::MpiWrapper& mpi);
        void RetrieveHamiltonian(utilities::MpiWrapper& mpi);
        void ClearHamiltonian();
        void EigensystemToFile(const bool writeEigenvalues, const bool writeEigenvectors, utilities::MpiWrapper& mpi);    
        void EigensystemFromFile(const bool readEigenvalues, const bool readEigenvectors, utilities::MpiWrapper& mpi);
    };
}   //  End diagonalization namespace
#endif
