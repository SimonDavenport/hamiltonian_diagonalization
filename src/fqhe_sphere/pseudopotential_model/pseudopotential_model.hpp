////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store a single Landau level FQHE 
//!     Haldane pseduopotential model for in the sphere geometry. 
//!
//!     The class allows for either dense or sparse matrix storage, and within 
//!     that it allows for a mapped sparse matrix or compressed column format.
//!     Any format can be selected to populate the matrix (though the mapped
//!     sparse matrix format is typically faster - and the dense format will
//!     use a lot more memory). The storage format will automatically
//!     be converted to CCS when any sparse diagonalization routine is called.  
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

#ifndef _PSEUDOPOTENTIAL_MODEL_HPP_INCLUDED_
#define _PSEUDOPOTENTIAL_MODEL_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "pseudopotential_model_base.hpp"
#include "../../hamiltonians/spinless_fermion_hamiltonian.hpp"
#include "term_tables.hpp"
#include "angular_momentum_constraint.hpp"
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#include <bitset>
#endif

namespace diagonalization
{
    //!
    //! The SpherePseudopotentialModel class defines a single Landau level FQHE
    //! Haldane pseudopotential model the sphere geometry
    //!
    class SpherePseudopotentialModel : 
    public SpherePseudopotentialModelBase<SpinlessFermionHamiltonian<double> >
    {
        private:
        QuadraticTermTables m_quadraticTables;
                                            //!<  Class container for regular array term tables
        QuarticTermTables m_quarticTables;  //!<  Class container for regular array term tables

        public:
        SpherePseudopotentialModel(const iSize_t nbrParticles, const iSize_t nbrOrbitals,
                                   const std::vector<double>& pseudopotentials);
        SpherePseudopotentialModel(boost::program_options::variables_map* optionList,
                                   utilities::MpiWrapper& mpi);
        ~SpherePseudopotentialModel();
        void BuildTermTables(const utilities::MpiWrapper& mpi);
        void TermsToFile(const io::fileFormat_t format, utilities::MpiWrapper& mpi) const;
        void TermsFromFile(const io::fileFormat_t format, utilities::MpiWrapper& mpi);
        void SetOccupationEnergies(std::vector<double>& energyLevels, const utilities::MpiWrapper& mpi);
        void BuildFockBasis(utilities::MpiWrapper& mpi);
        void SetFockBasis(fock_t* buffer, const fock_t dim);
        void BuildHamiltonian(utilities::MpiWrapper& mpi);
        void Diagonalize(utilities::MpiWrapper& mpi);
    };
}   //  End namespace diagonalization
#endif
