////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 13/02/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store a parameters defining a
//!     FQHE Haldane pseduopotential Hamiltonian in the sphere geometry
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

#ifndef _PSEUDOPOTENTIAL_HAMILTONIAN_DATA_HPP_INCLUDED_
#define _PSEUDOPOTENTIAL_HAMILTONIAN_DATA_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../utilities/general/cout_tools.hpp" 
                                        //  For cout manipulation
#include "../../utilities/wrappers/mpi_wrapper.hpp"   
                                        //  For MPI functionality
#include "../../utilities/mathematics/clebsch_gordan_coefficients.hpp"
                                        //  For Clebsch Gordan coefficients
#include <vector>
#include <boost/program_options.hpp>

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

namespace diagonalization
{

///////     ENUM AND STATIC DECLARATIONS     ///////////////////////////////////

enum io_t {_IN_,_OUT_};
//!<    Define a type for specifying construction of input or output files

enum diagonalizationMethod_t {_FULL_=0,_LANCZOS_=1};
//!<    Define a list of diagonalization methods:
//!     -FULL requires dense matrix storage and uses a QR decomposition
//!     -LANCZOS uses the ARPACK library and a compatible sparse matrix storage scheme

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief A data structure to contain all of the principle model parameters
//!
////////////////////////////////////////////////////////////////////////////////

struct PseudopotentialHamiltonianData
{
    iSize_t m_nbrOrbitals;              //!<    Number of angular momentum orbitals
    iSize_t m_nbrParticles;             //!<    Number of particles in the system
    mState_t m_maxLz;                   //!<    Twice the minimum Lz value used in orbital labelling the LLL
    mState_t m_maxLz2;                  //!<    Twice the minimum Lz value used in orbital labelling the 2nd LL
    mState_t m_totalLz;                 //!<    Total orbital angular momentum sector
    iSize_t m_nbrLevels;                //!<    Number of Landau levels
    std::vector<double> m_twoBodyPseudopotentials; 
                                        //!<    List of two-body pseudopotential coefficients
    std::vector<double> m_twoBodyPseudopotentials2LL; 
                                        //!<    List of 2LL two-body pseudopotential coefficients                                    
    std::string m_outPath;              //!<    Path to output data directory
	std::string m_inPath;               //!<    Path to input data directory
	std::string m_outFileName;          //!<    Name of file where output data are stored
	iSize_t m_nbrEigenvaluesToFind;    
	                                    //!<    A number to specify how many of the 
	                                    //!     lowest eigenvalues to find
	bool m_blockDiagonalize;            //!<    Set to true to diagonalize the Hamiltonian in 
                                        //!<    separate block-diagonal total angular momentum sectors
    diagonalizationMethod_t m_method;   //!<    Store the diagonalization method
    std::string m_initialVectorFile;    //!<    Name of ARPACK initial vector file
    std::string m_finalVectorFile;      //!<    Name of ARPACK final vector file  
    bool m_lookupTablesBuilt;           //!<    Set to true once the look-up tables are built
    bool m_fockBasisBuilt;              //!<    Set to true once the Fock basis has been constructed
    bool m_hamiltonianBuilt;            //!<    Set to true once the Hamiltonian is built
    bool m_hamiltonianDiagonalized;     //!<    Set to true once Hamiltonian is diagonalized

    //  Default Constructor
    PseudopotentialHamiltonianData();

    //  Constructor from minimal data
    PseudopotentialHamiltonianData(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals,const std::vector<double>& pseudopotentials);
    
    PseudopotentialHamiltonianData(const iSize_t nbrParticles,
        const iSize_t nbrOrbitals,const std::vector<double>& pseudopotentials,
        const std::vector<double>& pseudopotentials2LL);
    
    //  Destructor
    ~PseudopotentialHamiltonianData();

    //  Copy constructor
    PseudopotentialHamiltonianData(const PseudopotentialHamiltonianData& other);
    
    //  Constructor from command line arguments
    PseudopotentialHamiltonianData(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);
    
    //  Overload assignment operator
    PseudopotentialHamiltonianData& operator=(const PseudopotentialHamiltonianData& other);
    
    // MPI function to synchronise all variables with those on the master node
    void MpiSynchronize(const int nodeId,const utilities::MpiWrapper& mpi);

    //  Generate an output file name base, using struct parameters
    std::string MakeBaseFileName(const io_t io) const;
    
    //  Set the minimum value of Lz such that the orbitals are labelled by 
    //  -m_maxLz to m_maxLz in steps of 2
    void SetMaxLz(const iSize_t nbrOrbitals,const iSize_t nbrLevels); 
    
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End diagonalization namespace

#endif
