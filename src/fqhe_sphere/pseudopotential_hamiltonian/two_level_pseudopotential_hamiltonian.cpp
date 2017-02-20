////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 13/02/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store a FQHE Haldane pseudopotential
//!     Hamiltonian in the sphere geometry. Both 1 and 2LL systems
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "two_level_pseudopotential_hamiltonian.hpp"

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Constructor with just the number of particles, orbitals and
//! pseudopotentials specified
//!
////////////////////////////////////////////////////////////////////////////////

SphereTwoLevelPseudopotentialHamiltonian::SphereTwoLevelPseudopotentialHamiltonian(
    const iSize_t nbrParticles,                     //!<    Number of particles
    const iSize_t nbrOrbitals,                      //!<    Number of orbitals
    const std::vector<double>& pseudopotentials,    //!<    List of two-body pseudopotentials
    const std::vector<double>& pseudopotentials2LL) //!<    List of 2nd LL two-body pseudopotentials
{
    this->m_params = PseudopotentialHamiltonianData(nbrParticles,nbrOrbitals,pseudopotentials,pseudopotentials2LL);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Constructor from command line arguments
//!
////////////////////////////////////////////////////////////////////////////////

SphereTwoLevelPseudopotentialHamiltonian::SphereTwoLevelPseudopotentialHamiltonian(
    boost::program_options::variables_map* optionList,
                                //!<    Parsed command line argument list
    utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
:
    SpherePseudopotentialHamiltonianBase<TwoLevelSpinlessFermionHamiltonian<double> >(optionList,mpi)
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate look-up tables for momentum conservation and calculate
//! Vkkkk matrix elements
//! 
////////////////////////////////////////////////////////////////////////////////

void SphereTwoLevelPseudopotentialHamiltonian::BuildLookUpTables(
    const utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
{
    //////  POPULATE LOOK-UP TABLES FOR Vkkkk       ////////////////////////////

    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{ 
        utilities::cout.MainOutput()<<"\n\t============ BUILDING MATRIX ELEMENT LOOK-UP TABLES ============"<<std::endl;
    }
    
    //  Set # LLL orbitals
        
    iSize_t nbrOrbitalsLLL = (m_params.m_nbrOrbitals-2)/2;
    
    //  Allocate memory to store the k-value Vkkkk value tables and  
    //  populate the tables

    m_quarticTables.Initialize(nbrOrbitalsLLL,mpi);

    this->BaseGenerate4KTable(m_quarticTables.GetKTable(),0,mpi);

    this->BaseGetCoefficientsFromPseudopotentials(m_quarticTables.GetVTable(),0,mpi);
    
    #if _DEBUG_
    PAR_PRINT("LLL",m_quarticTables.GetVTable()->data(),m_quarticTables.GetVTable()->size());
    PAR_PAUSE();
    #endif
    
    m_quarticTables2LL.Initialize(nbrOrbitalsLLL+2,mpi);
    
    this->BaseGenerate4KTable(m_quarticTables2LL.GetKTable(),1,mpi);

    this->BaseGetCoefficientsFromPseudopotentials(m_quarticTables2LL.GetVTable(),1,mpi);

    #if _DEBUG_
    PAR_PRINT("2LL",m_quarticTables2LL.GetVTable()->data(),m_quarticTables2LL.GetVTable()->size());
    PAR_PAUSE();
    #endif

    //  Synchronize with the master node
       
    if(0 == mpi.m_id)    // FOR THE MASTER NODE
    { 
        utilities::cout.AdditionalInfo()<<"\n\t- MPI SYNC HAMILTONIAN COEFFICIENT TABLES"<<std::endl;
    }

    m_quarticTables.MpiSynchronize(0,mpi);
    m_quarticTables2LL.MpiSynchronize(0,mpi);

    m_params.m_lookupTablesBuilt = true;

    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Set optional CdC terms
//! 
////////////////////////////////////////////////////////////////////////////////

void SphereTwoLevelPseudopotentialHamiltonian::SetOccupationEnergies(
    double* energyLevels,               //!<    Single particle energy levels
    const iSize_t dim,                  //!<    Length of array given
    const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{
    if(dim != m_params.m_nbrOrbitals)
    {
        std::cerr<<"\n\tERROR in SetOccupationEnergies: Must specify same number of energy levels as orbitals"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    //  Set # LLL orbitals
        
    iSize_t nbrOrbitalsLLL = (m_params.m_nbrOrbitals-2)/2;
    
    m_quadraticTables.Initialize(nbrOrbitalsLLL,mpi);
    m_quadraticTables2LL.Initialize(nbrOrbitalsLLL+2,mpi);
    
    memcpy(m_quadraticTables.GetVTable()->data(),energyLevels,sizeof(double)*nbrOrbitalsLLL);
    memcpy(m_quadraticTables2LL.GetVTable()->data(),energyLevels+nbrOrbitalsLLL,sizeof(double)*(nbrOrbitalsLLL+2));
    
    //PRINT("energyLevels",m_quadraticTables.GetVTable()->data(),m_quadraticTables.GetVTable()->size());
    //PRINT("energyLevels",m_quadraticTables2LL.GetVTable()->data(),m_quadraticTables2LL.GetVTable()->size());
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate a Fock basis 
//! 
////////////////////////////////////////////////////////////////////////////////

void SphereTwoLevelPseudopotentialHamiltonian::BuildFockBasis(
    utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
{
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {   
        if(m_params.m_blockDiagonalize)
        {
            utilities::cout.MainOutput()<<"\n\t============ BUILDING FOCK BASIS (lz_tot = "<<m_params.m_totalLz<<") ============ "<<std::endl;
        }
        else
        {
            utilities::cout.MainOutput()<<"\n\t============ BUILDING FULL FOCK BASIS ============ "<<std::endl;
        }
    }
    
    //  The function "GenerateFockSpace" requires a single argument:
    //  a function that takes a Fock basis state as an argument and returns 
    //  a bool saying whether that state is in the specified sector or not. 

    bool isZeroDimensional = false;
    
    iSize_t nbrOrbitalsLLL = (m_params.m_nbrOrbitals-2)/2;

    if(m_params.m_blockDiagonalize)
    {
        TwoLevelFermionFockBasis basis;
    
        //  To make the function we're going to use std::bind to fix parameters 
        //  in the function called TestLinearMomentumSector defined above, such 
        //  that only the first parameter (the Fock state) can be 
        //  specified and all other parameters are fixed.

        isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles,nbrOrbitalsLLL,nbrOrbitalsLLL+2,std::bind(&TestAngularMomentumSector2LL,std::placeholders::_1,m_params.m_totalLz,m_params.m_maxLz2),mpi);

        m_hamiltonian.SetFockBasis(basis);
    }
    else
    {
        TwoLevelFermionFockBasis basis;
    
        //  If we're not block diagonalizing then we pass a dummy function
        //  to m_hamiltonian (the function is declared as a "lambda function")
        //  This function will always return true.
        
        isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles,nbrOrbitalsLLL,nbrOrbitalsLLL+2,[](const TwoLevelFockState& state){return true;},mpi);
        
        m_hamiltonian.SetFockBasis(basis);
    }

    if(isZeroDimensional)
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        { 
            utilities::cout.SecondaryOutput()<<"\n\t============ FOCK BASIS IS ZERO DIMENSIONAL ============ "<<std::endl;
        }
        
        m_params.m_fockBasisBuilt = false;
    
        return;
    }
    
    m_params.m_fockBasisBuilt = true;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Set Fock basis from external array
//!
////////////////////////////////////////////////////////////////////////////////

void SphereTwoLevelPseudopotentialHamiltonian::SetFockBasis(
    fock_t* buffer1,        //!<    Array to set lower LLs from
    fock_t* buffer2,        //!<    Array to set upper LLs from
    const fock_t dim)       //!<    Dimension of the array
{
    TwoLevelFermionFockBasis basis;
    
    basis.SetFockBasis(buffer1,buffer2,dim,0);

    m_hamiltonian.SetFockBasis(basis);
    
    m_params.m_fockBasisBuilt = true;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate a matrix representation of our Hamiltonian in the 2D 
//! momentum space basis.
//!
//! When called multiple times, this function will not construct a new 
//! unless the previous Hamiltonian has been cleared with the 
//! ClearHamiltonian function.
//! 
////////////////////////////////////////////////////////////////////////////////

void SphereTwoLevelPseudopotentialHamiltonian::BuildHamiltonian(
    utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper
{
    if(m_params.m_lookupTablesBuilt && m_params.m_fockBasisBuilt && !m_params.m_hamiltonianBuilt)
    {
        //////  GENERATE THE HAMILTONIAN        //////
        
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {   
            utilities::cout.MainOutput()<<"\n\t============ BUILDING HAMILTONIAN ============ "<<std::endl;
        }

        m_hamiltonian.Initialize(m_params.m_nbrParticles,
        m_params.m_nbrOrbitals,utilities::_SPARSE_MAPPED_,mpi);

        MPI_Barrier(mpi.m_comm);

        //////      Build the Hamiltonian

        if(m_quadraticTables.GetDimension()>0)
        {
            m_hamiltonian.Add_CdC_Terms(&m_quadraticTables,0,0,0,m_hamiltonian.m_data.m_fockSpaceDim,mpi);
        }
        
        if(m_quadraticTables2LL.GetDimension()>0)
        {
            m_hamiltonian.Add_CdC_Terms(&m_quadraticTables2LL,1,1,0,m_hamiltonian.m_data.m_fockSpaceDim,mpi);
        }
        
        //m_hamiltonian.PrintHamiltonian(0);
        
        m_hamiltonian.Add_CdCdCC_Terms(&m_quarticTables,0,0,0,0,0,m_hamiltonian.m_data.m_fockSpaceDim,mpi);
        
        m_hamiltonian.Add_CdCdCC_Terms(&m_quarticTables2LL,1,1,1,1,0,m_hamiltonian.m_data.m_fockSpaceDim,mpi);
                
        m_hamiltonian.PrintMemoryAllocation(mpi);

        m_hamiltonian.PrintHamiltonian(0,mpi);

        m_params.m_hamiltonianBuilt = true;

        //std::stringstream fileName;
        //fileName.str("");
        //fileName<<"temp_"<<mpi.m_id<<"_hamiltonian.dat";
        //m_hamiltonian.HamiltonianToFile(fileName.str(),mpi);
    }
    else
    {
        m_params.m_hamiltonianBuilt = false;
    }
}
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Public interface to diagonalize the Hamiltonian
//!
//! TODO: update and complete this implementation to apply correctly to 
//! the two-level system
//!
////////////////////////////////////////////////////////////////////////////////
    
void SphereTwoLevelPseudopotentialHamiltonian::Diagonalize(
    utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
{
    this->BaseDiagonalize(&m_quadraticTables,&m_quarticTables,mpi);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace diagonalization 

