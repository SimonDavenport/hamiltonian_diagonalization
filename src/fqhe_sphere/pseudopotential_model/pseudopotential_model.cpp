////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store a single Landau level FQHE 
//!     Haldane pseduopotential model for in the sphere geometry. 
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "pseudopotential_model.hpp"

namespace diagonalization
{
    //!
    //! Constructor with just the number of particles, orbitals and
    //! pseudopotentials specified
    //!
    SpherePseudopotentialModel::SpherePseudopotentialModel(
        const iSize_t nbrParticles,                     //!<    Number of particles
        const iSize_t nbrOrbitals,                      //!<    Number of orbitals
        const std::vector<double>& pseudopotentials)    //!<    List of two-body pseudopotentials
    {
        this->m_params = PseudopotentialModelData(nbrParticles, nbrOrbitals, pseudopotentials);
    }

    //!
    //! Constructor from command line arguments
    //!
    SpherePseudopotentialModel::SpherePseudopotentialModel(
        boost::program_options::variables_map* optionList,
                                    //!<    Parsed command line argument list
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    :
        SpherePseudopotentialModelBase<SpinlessFermionHamiltonian<double> >(optionList, mpi)
    {}

    //!
    //! Destructor
    //!
    SpherePseudopotentialModel::~SpherePseudopotentialModel()
    {}

    //!
    //! Generate term tables for momentum conservation and calculate
    //! Vkkkk matrix elements
    //! 
    void SpherePseudopotentialModel::BuildTermTables(
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        //////  POPULATE TABLES FOR Vkkkk       ////////////////////////////////
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t============ BUILDING TERM LOOK-UP TABLES ============"<<std::endl;
        }
        //  Allocate memory to store the k-value Vkkkk value tables and populate 
        //  the tables
        m_quarticTables.Initialize(m_params.m_nbrOrbitals, mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING MOMENTUM CONSERVING k1,k2,k3,k4 TABLE"<<std::endl;
        }
        this->BaseGenerate4KTable(m_quarticTables.GetKTable(), 0, mpi);
	    MPI_Barrier(mpi.m_comm);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING 2-BODY PSEUDOPOTENTIAL COEFFICIENTS"<<std::endl;
        }
        this->BaseGetCoefficientsFromPseudopotentials(m_quarticTables.GetVTable(), 0, mpi);
	    MPI_Barrier(mpi.m_comm);
        //  Synchronize with the master node
        if(0 == mpi.m_id)    // FOR THE MASTER NODE
        { 
            utilities::cout.AdditionalInfo()<<"\n\t- MPI SYNC TERM TABLES"<<std::endl;
        }
        m_quarticTables.MpiSynchronize(0, mpi);
        m_params.m_lookupTablesBuilt = true;
        return;
    }
    
    //!
    //! Store term tables in a file
    //!
    SpherePseudopotentialModel::TermTablesToFile(
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t============ WRITING TERM LOOK-UP TABLES TO FILE ============"<<std::endl;
        }
        m_quadraticTables.ToFile();
        m_quarticTables.ToFile();
        return;
    }
    
    //!
    //! Retrieve term tables from a file
    //!
    SpherePseudopotentialModel::TermTablesFromFile(
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t============ RETRIEVING TERM LOOK-UP TABLES FROM FILE ============"<<std::endl;
        }
        m_quadraticTables.FromFile();
        m_quarticTables.FromFile();
        return;
    }
    
    //!
    //! Set optional CdC terms
    //! 
    void SpherePseudopotentialModel::SetOccupationEnergies(
        double* energyLevels,               //!<    Single particle energy levels
        const iSize_t dim,                  //!<    Length of array given
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(dim != m_params.m_nbrOrbitals)
        {
            std::cerr<<"\n\tERROR in SetOccupationEnergies: Must specify same number of energy levels as orbitals"<<std::endl;
            exit(EXIT_FAILURE);
        }
        m_quadraticTables.Initialize(m_params.m_nbrOrbitals, mpi);
        memcpy(m_quadraticTables.GetVTable()->data(), energyLevels, sizeof(double)*dim);
    }

    //!
    //! Generate a Fock basis 
    //! 
    void SpherePseudopotentialModel::BuildFockBasis(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
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
        if(m_params.m_blockDiagonalize)
        {
            FermionFockBasis basis;
            //  To make the function we're going to use std::bind to fix parameters 
            //  in the function called TestLinearMomentumSector defined above, such 
            //  that only the first parameter (the Fock state) can be 
            //  specified and all other parameters are fixed.
            isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, m_params.m_nbrOrbitals, 
                                                        std::bind(&TestAngularMomentumSector, std::placeholders::_1,
                                                        m_params.m_totalLz, m_params.m_maxLz), mpi);
            m_hamiltonian.SetFockBasis(basis);
        }
        else
        {
            FermionFockBasis basis;
            //  If we're not block diagonalizing then we pass a dummy function
            //  to m_hamiltonian (the function is declared as a "lambda function")
            //  This function will always return true.
            isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, m_params.m_nbrOrbitals, 
                                                        [](const fock_t state){return true;}, mpi);
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

    //!
    //! Set Fock basis from external array
    //!
    void SpherePseudopotentialModel::SetFockBasis(
        fock_t* buffer,         //!<    Array to set from
        const fock_t dim)       //!<    Dimension of the array
    {
        FermionFockBasis basis;
        basis.SetFockBasis(buffer, dim, 0);
        m_hamiltonian.SetFockBasis(basis);
        m_params.m_fockBasisBuilt = true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a matrix representation of our Hamiltonian in the 2D momentum
    //! space basis.
    //!
    //! When called multiple times, this function will not construct a new Hamiltonian
    //! unless the previous Hamiltonian has been cleared with the ClearHamiltonian
    //! function.
    ////////////////////////////////////////////////////////////////////////////////
    void SpherePseudopotentialModel::BuildHamiltonian(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(m_params.m_lookupTablesBuilt && m_params.m_fockBasisBuilt && !m_params.m_hamiltonianBuilt)
        {
            //////  GENERATE THE HAMILTONIAN        //////
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {   
                utilities::cout.MainOutput()<<"\n\t============ BUILDING HAMILTONIAN ============ "<<std::endl;
            }
            m_hamiltonian.Initialize(m_params.m_nbrParticles, m_params.m_nbrOrbitals, utilities::_SPARSE_MAPPED_, mpi);
            MPI_Barrier(mpi.m_comm);
            //////      Build the Hamiltonian
            if(m_quadraticTables.GetDimension()>0)
            {
                m_hamiltonian.Add_CdC_Terms(&m_quadraticTables, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
            }
            m_hamiltonian.Add_CdCdCC_Terms(&m_quarticTables, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
            m_hamiltonian.PrintMemoryAllocation(mpi);
            m_hamiltonian.PrintHamiltonian(0, mpi);
            m_params.m_hamiltonianBuilt = true;
        }
        else
        {
            m_params.m_hamiltonianBuilt = false;
        }
    }

    //!
    //! Public interface to diagonalize the Hamiltonian
    //! 
    void SpherePseudopotentialModel::Diagonalize(
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    {
        this->BaseDiagonalize(&m_quadraticTables, &m_quarticTables, mpi);
    }
}   //  End namespace diagonalization