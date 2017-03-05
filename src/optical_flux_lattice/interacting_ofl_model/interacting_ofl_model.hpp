////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store the optical flux lattice interacting
//!     model. See e.g. PRL 109, 265301(2012) 	   
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

#ifndef _INTERACTING_OFL_MODEL_HPP_INCLUDED_
#define _INTERACTING_OFL_MODEL_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include <functional>
#include <string.h>
#include "../program_options/interacting_ofl_model_options.hpp"
#include "../program_options/noninteracting_ofl_model_options.hpp"
#include "../noninteracting_ofl_model/noninteracting_ofl_model.hpp"
#include "../noninteracting_ofl_model/noninteracting_ofl_model_grid.hpp"
#include "interacting_ofl_model_data.hpp"
#include "../../hamiltonians/spinless_fermion_hamiltonian.hpp"
#include "../../observables/fermion_observables.hpp"
#include "term_tables.hpp"
#include "term_hash_tables.hpp"
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
#include "../../utilities/wrappers/program_options_wrapper.hpp"
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#include <bitset>
#endif

namespace diagonalization
{
    //!
    //! The InteractingOflModel class contains a matrix 
    //! representation of the interacting optical flux lattice model
    //!
    class InteractingOflModel
    {
        friend class ObservablesManager;
        private:
        SpinlessFermionHamiltonian<dcmplx> m_hamiltonian;   
                                            //!<  An object to contain the Hamiltonian itself
        InteractingOflModelData m_params;   //!<  Class internal implementation of the     
                                            //!   InteractingOflModelData struct          
        NonInteractingOflModelData m_oflParameters;
                                            //!<  Instance of the non-interacting 
                                            //!   ofl parameters class
        QuadraticTermTables m_quadraticTables;
                                            //!<  Class container for regular array look-up tables
                                            //!   used in the case where we have full 2D momentum conservation 
        QuarticTermTables m_quarticTables;        
                                            //!<  Class container for regular array look-up tables
                                            //!   used in the case where we have full 2D momentum conservation 
        QuadraticTermHashTables m_quadraticHashTables;
                                            //!<  Class container for multi-hash maps:
                                            //!   used when we have only ky-momentum conservation
        QuarticTermHashTables m_quarticHashTables;
                                            //!<  Class container for multi-hash maps:
                                            //!   used when we have only ky-momentum conservation
        std::vector<dcmplx> m_basisChangeMatrix;        
                                            //!<    Change of basis matrix;
        std::vector<double> m_magnetization;//!<  A list of the magnetization values
        
        void Generate4KTable(std::vector<kState_t>* kTable,std::vector<int>* gxTable,std::vector<int>* gyTable, 
                             const utilities::MpiWrapper& mpi);    
        void GenerateWannier2KTable(utilities::MultiHashMultiMap<kState_t>* kHashTable, const utilities::MpiWrapper& mpi);
        void GenerateWannier4KTable(utilities::MultiHashMultiMap<kState_t>* kHashTable, 
                                    const utilities::MpiWrapper& mpi);
        kState_t FindK1(const kState_t k2, const kState_t k3, const kState_t k4, int& gTot1, int& gTot2) const;
        kState_t FindK1(const kState_t k2y, const kState_t k3y, const kState_t k4y, int& gTot2) const;
        void GenerateQuarticTerms(std::vector<dcmplx>* table, const iSize_t dimension, NonInteractingOflModel* blochTable, 
                                  std::vector<int>* gxTable, std::vector<int>* gyTable, utilities::MpiWrapper& mpi) const;
        void ConvertTableFormat(const utilities::MpiWrapper& mpi);
        void ChangeToWannierBasis(NonInteractingOflModel* blochTable, const double matrixElementTol, 
                                  utilities::MpiWrapper& mpi);
	    void CalcualteSigmaZMap(NonInteractingOflModel* blochTable, const utilities::MpiWrapper& mpi);
        bool CalculateOccupations(utilities::MpiWrapper& mpi) const;
        bool CalculateSusceptibility(utilities::MpiWrapper& mpi) const;
        bool PlotHamiltonian(utilities::MpiWrapper& mpi);
        bool ListMostProbable(const iSize_t nbrStates,utilities::MpiWrapper& mpi) const;
        bool CalculateDensityDensityFunction(utilities::MpiWrapper& mpi) const;
        bool CalculateParticipationRatio(utilities::MpiWrapper& mpi) const;
        bool CalculateTranslationalDensityDensityFunction(utilities::MpiWrapper& mpi) const;
        bool CalculateRotationalDensityDensityFunction(utilities::MpiWrapper& mpi) const;
        public:
        InteractingOflModel();  
        InteractingOflModel(InteractingOflModelData params);
        InteractingOflModel(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        ~InteractingOflModel();
        void UpdateSqlStatus(const utilities::MpiWrapper& mpi);
        void BuildTermTables(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void TermsToFile(const io::fileFormat_t format, utilities::MpiWrapper& mpi);
	    void TermsFromFile(const io::fileFormat_t format, utilities::MpiWrapper& mpi);
        void SetSector(const kState_t kxTot,const kState_t kyTot);
        void BuildFockBasis(utilities::MpiWrapper& mpi);
        void BuildHamiltonian(utilities::MpiWrapper& mpi);
        void Diagonalize(utilities::MpiWrapper& mpi);
        void HamiltonianToFile(const io::fileFormat_t format, utilities::MpiWrapper& mpi);
        void HamiltonianFromFile(const io::fileFormat_t format, utilities::MpiWrapper& mpi);
        void ClearHamiltonian();
        void EigensystemToFile(const bool writeEigenvalues, const bool writeEigenvectors, 
                               const io::fileFormat_t format, utilities::MpiWrapper& mpi);    
        void EigensystemFromFile(const bool readEigenvalues, const bool readEigenvectors, 
                                 const io::fileFormat_t format, utilities::MpiWrapper& mpi);
    };
}   //  End diagonalization namespace
#endif
