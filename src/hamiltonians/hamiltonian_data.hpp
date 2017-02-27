////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains an underlying data structure to store variables 
//!     in the Hamiltonian class
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

#ifndef _FERMION_HAMILTONIAN_DATA_HPP_INCLUDED_
#define _FERMION_HAMILTONIAN_DATA_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../utilities/data_structures/sparse_matrix.hpp"
#include "../utilities/general/orbital_and_state_defs.hpp"
#include "../program_options/arpack_options.hpp"
#include "../utilities/general/dcmplx_type_def.hpp" 
#include "../utilities/wrappers/mpi_wrapper.hpp"
#include "../utilities/wrappers/arpack_wrapper.hpp" 
#include "../utilities/general/cout_tools.hpp"
#include "../utilities/mathematics/binary_number_tools.hpp"
#include <boost/program_options.hpp>
#include <vector>
#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    //////      Declare some utility functions      ////////////////////////////////
    void ValueToFile(std::ofstream& f_out, double* value);
    void ValueToFile(std::ofstream& f_out, dcmplx* value);
    void ValueFromFile(std::ifstream& f_in, double* value);
    void ValueFromFile(std::ifstream& f_in, dcmplx* value);
    ///////     ENUM AND STATIC DECLARATIONS     ///////////////////////////////////
    enum parallelFlag_t {_PARALLEL_,_SERIAL_};
    static const fock_t _MEMORY_QUERY_LIMIT_ = 4294967296;
    static constexpr double _SPARSE_ZERO_TOL_ = 0.00000000000001;
    //!
    //! A data structure to contain all of the principle model parameters
    //!
    struct HamiltonianData
    {
        iSize_t m_nbrParticles;             //!<    Number of particles described
        iSize_t m_nbrOrbitals;              //!<    Number of orbitals  
        fock_t m_highestOrbital;            //!<    Highest allowed k-state in its 
                                            //!     Fock-state representation
        utilities::storageMethod_t m_storageMethod;	
                                            //!<	Flag specifying if we store the matrix 
                                            //!     in sparse or dense format
        bool m_matrixAllocated;             //!<    Flag to say whether the matrix is allocated or not
        bool m_eigenvaluesCalculated;       //!<    Flag to say whether the eigenvalues 
                                            //!     are calculated or not
        bool m_eigenvectorsCalculated;      //!<    Flag to say whether the eigenvalues 
                                            //!     are calculated or not
        fock_t m_fockSpaceDim;              //!<    Total dimension of the Fock space
        fock_t m_nodeDim;                   //!<    Total dimension of the Fock space
                                            //!     stored on a given node
        iSize_t m_nbrEigenvalues;           //!<    Number of eigenvalues to calculate
        HamiltonianData();
        HamiltonianData(const iSize_t nbrParticles,
                        const iSize_t nbrOrbitals, const utilities::storageMethod_t storageMethod);
        HamiltonianData(const HamiltonianData& other);
        HamiltonianData& operator=(const HamiltonianData& other);
        ~HamiltonianData();
        void Clear();
        void SetDataDimensions(const fock_t fockSpaceDimension,
                               const iSize_t nbrNodes, utilities::MpiWrapper& mpi);
    };
}   //  End namespace diagonalization 
#endif
