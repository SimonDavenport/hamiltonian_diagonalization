////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store a FQHE Haldane pseudopotential
//!     model in the sphere geometry. 
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
#include "two_level_pseudopotential_model.hpp"

namespace diagonalization
{
    //!
    //! Constructor with just the number of particles, orbitals and
    //! pseudopotentials specified
    //!
    SphereTwoLevelPseudopotentialModel::SphereTwoLevelPseudopotentialModel(
        const iSize_t nbrParticles,                     //!<    Number of particles
        const iSize_t nbrOrbitals,                      //!<    Number of orbitals
        const std::vector<double>& pseudopotentials,    //!<    List of two-body pseudopotentials
        const std::vector<double>& pseudopotentials2LL) //!<    List of 2nd LL two-body pseudopotentials
    {
        this->m_params = PseudopotentialModelData(nbrParticles, nbrOrbitals,
                                                  pseudopotentials, pseudopotentials2LL);
    }

    //!
    //! Constructor from command line arguments
    //!
    SphereTwoLevelPseudopotentialModel::SphereTwoLevelPseudopotentialModel(
        boost::program_options::variables_map* optionList,
                                    //!<    Parsed command line argument list
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    :
        SpherePseudopotentialModelBase<TwoLevelSpinlessFermionHamiltonian<double> >(optionList, mpi)
    {}

    //!
    //! Generate look-up tables for momentum conservation and calculate
    //! Vkkkk matrix elements
    //!

    void SphereTwoLevelPseudopotentialModel::BuildTermTables(
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
        m_quarticTables.Initialize(nbrOrbitalsLLL, mpi);
        this->BaseGenerate4KTable(m_quarticTables.GetKTable(), 0, mpi);
        this->BaseTermsFromPseudopotentials(m_quarticTables.GetVTable(), 0, mpi);
        m_quarticTables2LL.Initialize(nbrOrbitalsLLL+2, mpi);
        this->BaseGenerate4KTable(m_quarticTables2LL.GetKTable(), 1, mpi);
        this->BaseTermsFromPseudopotentials(m_quarticTables2LL.GetVTable(), 1, mpi);
        if(m_params.m_background != 0.0)
	    {
	        std::vector<double> backgroundPerParticle(m_params.m_nbrOrbitals, 
	                                                  -this->m_params.m_background / m_params.m_nbrParticles);
	        this->SetOccupationEnergies(backgroundPerParticle, mpi);
	    }
        //  Synchronize with the master node
        if(0 == mpi.m_id)    // FOR THE MASTER NODE
        { 
            utilities::cout.AdditionalInfo()<<"\n\t- MPI SYNC HAMILTONIAN COEFFICIENT TABLES"<<std::endl;
        }
        m_quarticTables.MpiSynchronize(0, mpi);
        m_quarticTables2LL.MpiSynchronize(0, mpi);
        if(_HASH_ == m_params.m_setTableFormat)
        {
            this->ConvertTableFormat(mpi);
        }
        m_params.m_termTablesBuilt = true;
        return;
    }
    
    //!
    //! This function converts the table format for storing Vkkkk, Ekk and 
    //! the momentum conserving tables from a regular array format to a multi
    //! key hash table format, if not done so already
    //!
    void SphereTwoLevelPseudopotentialModel::ConvertTableFormat(
        const utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class 
    {
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t- CONVERTING LOOK-UP ARRAYS TO HASH TABLES"<<std::endl;
            }
            m_quadraticHashTables.Clear();
            m_quarticHashTables.Clear();
            m_quadraticHashTables2LL.Clear();
            m_quarticHashTables2LL.Clear();
            m_quadraticHashTables.Initialize(m_params.m_nbrOrbitals);
            m_quarticHashTables.Initialize(m_params.m_nbrOrbitals);
            m_quadraticHashTables2LL.Initialize(m_params.m_nbrOrbitals);
            m_quarticHashTables2LL.Initialize(m_params.m_nbrOrbitals);
            m_quadraticHashTables.SetFromArray(&m_quadraticTables, m_params.m_nbrOrbitals);
            m_quarticHashTables.SetFromArray(&m_quarticTables, m_params.m_nbrOrbitals);
            m_quadraticHashTables2LL.SetFromArray(&m_quadraticTables2LL, m_params.m_nbrOrbitals);
            m_quarticHashTables2LL.SetFromArray(&m_quarticTables2LL, m_params.m_nbrOrbitals);
            m_quarticTables.Clear();
            m_quadraticTables.Clear();
            m_quarticTables2LL.Clear();
            m_quadraticTables2LL.Clear();
            m_params.m_tableFormat = _HASH_;
            m_quadraticHashTables.MpiSynchronize(0, mpi);
            m_quarticHashTables.MpiSynchronize(0, mpi);
            m_quadraticHashTables2LL.MpiSynchronize(0, mpi);
            m_quarticHashTables2LL.MpiSynchronize(0, mpi);
        }
        return;
    }

    //!
    //! Set optional CdC terms
    //! 
    void SphereTwoLevelPseudopotentialModel::SetOccupationEnergies(
        std::vector<double>& energyLevels,  //!<    Single particle energy levels
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(energyLevels.size() != m_params.m_nbrOrbitals)
        {
            std::cerr<<"\n\tERROR in SetOccupationEnergies: Must specify same number of energy levels as orbitals"<<std::endl;
            exit(EXIT_FAILURE);
        }
        //  Set # LLL orbitals  
        iSize_t nbrOrbitalsLLL = (m_params.m_nbrOrbitals-2)/2;
        m_quadraticTables.Initialize(nbrOrbitalsLLL, mpi);
        m_quadraticTables2LL.Initialize(nbrOrbitalsLLL+2, mpi);
        memcpy(m_quadraticTables.GetVTable()->data(), energyLevels.data(), sizeof(double)*nbrOrbitalsLLL);
        memcpy(m_quadraticTables2LL.GetVTable()->data(), &energyLevels[0]+nbrOrbitalsLLL, sizeof(double)*(nbrOrbitalsLLL+2));
        if(_ARRAY_ != m_params.m_tableFormat)
        {
            m_quadraticHashTables.Initialize(m_params.m_nbrOrbitals);
            m_quadraticHashTables2LL.Initialize(m_params.m_nbrOrbitals);
            m_quadraticHashTables.SetFromArray(&m_quadraticTables, m_params.m_nbrOrbitals);
            m_quadraticHashTables2LL.SetFromArray(&m_quadraticTables2LL, m_params.m_nbrOrbitals);
            m_quadraticTables.Clear();
            m_quadraticTables2LL.Clear();
            m_params.m_tableFormat = _HASH_;
        }
    }

    //!
    //! Generate a Fock basis 
    //! 
    void SphereTwoLevelPseudopotentialModel::BuildFockBasis(
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
            isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, nbrOrbitalsLLL, nbrOrbitalsLLL+2,
                                                        std::bind(&TestAngularMomentumSector2LL, std::placeholders::_1, 
                                                        m_params.m_totalLz, m_params.m_maxLz2), mpi);
            m_hamiltonian.SetFockBasis(basis);
        }
        else
        {
            TwoLevelFermionFockBasis basis;
            //  If we're not block diagonalizing then we pass a dummy function
            //  to m_hamiltonian (the function is declared as a "lambda function")
            //  This function will always return true.
            isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, nbrOrbitalsLLL, nbrOrbitalsLLL+2,
                                                        [](const TwoLevelFockState& state){return true;}, mpi);
            m_hamiltonian.SetFockBasis(basis);
        }

        if(isZeroDimensional)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            { 
                utilities::cout.SecondaryOutput()<<"\n\tERROR: FOCK BASIS IS ZERO DIMENSIONAL ON ONE OR MORE NODES "<<std::endl;
            }
            m_params.m_fockBasisBuilt = false;
            return;
        }
        m_params.m_fockBasisBuilt = true;
    }

    //!
    //! Set Fock basis from external array
    //!
    void SphereTwoLevelPseudopotentialModel::SetFockBasis(
        fock_t* buffer1,        //!<    Array to set lower LLs from
        fock_t* buffer2,        //!<    Array to set upper LLs from
        const fock_t dim)       //!<    Dimension of the array
    {
        TwoLevelFermionFockBasis basis;
        basis.SetFockBasis(buffer1, buffer2, dim, 0);
        m_hamiltonian.SetFockBasis(basis);
        m_params.m_fockBasisBuilt = true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a matrix representation of our Hamiltonian in the 2D 
    //! momentum space basis.
    //!
    //! When called multiple times, this function will not construct a new 
    //! unless the previous Hamiltonian has been cleared with the 
    //! ClearHamiltonian function.
    ////////////////////////////////////////////////////////////////////////////////
    void SphereTwoLevelPseudopotentialModel::BuildHamiltonian(
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper
    {
        if(m_params.m_termTablesBuilt && m_params.m_fockBasisBuilt && !m_params.m_hamiltonianBuilt)
        {
            //////  GENERATE THE HAMILTONIAN        //////
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {   
                utilities::cout.MainOutput()<<"\n\t============ BUILDING HAMILTONIAN ============ "<<std::endl;
            }
            m_hamiltonian.Initialize(m_params.m_nbrParticles, m_params.m_nbrOrbitals, utilities::_SPARSE_MAPPED_, mpi);
            MPI_Barrier(mpi.m_comm);
            //////      Build the Hamiltonian
            if(_ARRAY_ == m_params.m_tableFormat)
            {
                if(m_quadraticTables.GetDimension())
                {
                    m_hamiltonian.Add_CdC_Terms(&m_quadraticTables, 0, 0, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                }
                if(m_quadraticTables2LL.GetDimension())
                {
                    m_hamiltonian.Add_CdC_Terms(&m_quadraticTables2LL, 0, 0, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                }
                m_hamiltonian.Add_CdCdCC_Terms(&m_quarticTables, 0, 0, 0, 0, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                m_hamiltonian.Add_CdCdCC_Terms(&m_quarticTables2LL, 1, 1, 1, 1, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi); 
            }
            else
            {
                if(m_quadraticHashTables.GetDimension())
                {
                    m_hamiltonian.Add_CdC_Terms(&m_quadraticHashTables, 0, 0, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                }
                if(m_quadraticHashTables2LL.GetDimension())
                {
                    m_hamiltonian.Add_CdC_Terms(&m_quadraticHashTables2LL, 0, 0, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                }
                m_hamiltonian.Add_CdCdCC_Terms(&m_quarticHashTables, 0, 0, 0, 0, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                m_hamiltonian.Add_CdCdCC_Terms(&m_quarticHashTables2LL, 1, 1, 1, 1, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi); 
            }   
            m_hamiltonian.PrintMemoryAllocation(mpi);
            //m_hamiltonian.PrintHamiltonian(0, mpi);
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
    void SphereTwoLevelPseudopotentialModel::Diagonalize(
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    {
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            this->BaseDiagonalize(&m_quadraticTables, &m_quarticTables, mpi);
        }
        else
        {
            this->BaseDiagonalize(&m_quadraticHashTables, &m_quarticHashTables, mpi);
        }
    }
    
    //!
    //! Store term tables in a file
    //!
    void SphereTwoLevelPseudopotentialModel::TermsToFile(
        const io::fileFormat_t format,      //!<    Format of file
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t============ WRITING TERM LOOK-UP TABLES TO FILE ============"<<std::endl;
        }
        std::stringstream fileNameQuadratic, filenameQuartic;
        fileNameQuadratic.str("");
        filenameQuartic.str("");
        fileNameQuadratic << m_params.m_outPath << "/quadratic_coefficient_table_LLL_L_" << m_params.m_nbrOrbitals << ".dat";
        filenameQuartic << m_params.m_outPath << "/quartic_coefficient_table_LLL_L_" << m_params.m_nbrOrbitals << ".dat";
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            m_quadraticTables.ToFile(fileNameQuadratic.str(), format, mpi);
            m_quarticTables.ToFile(filenameQuartic.str(), format, mpi);
        }
        else
        {
            m_quadraticHashTables.ToFile(fileNameQuadratic.str(), format, mpi);
            m_quarticHashTables.ToFile(filenameQuartic.str(), format, mpi);
        }
        // 2LL terms
        fileNameQuadratic.str("");
        filenameQuartic.str("");
        fileNameQuadratic << m_params.m_outPath << "/quadratic_coefficient_table_2LL_L_" << m_params.m_nbrOrbitals << ".dat";
        filenameQuartic << m_params.m_outPath << "/quartic_coefficient_table_2LL_L_" << m_params.m_nbrOrbitals << ".dat";
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            m_quadraticTables2LL.ToFile(fileNameQuadratic.str(), format, mpi);
            m_quarticTables2LL.ToFile(filenameQuartic.str(), format, mpi);
        }
        else
        {
            m_quadraticHashTables2LL.ToFile(fileNameQuadratic.str(), format, mpi);
            m_quarticHashTables2LL.ToFile(filenameQuartic.str(), format, mpi);
        }
        return;
    }
    
    //!
    //! Retrieve term tables from a file
    //!
    void SphereTwoLevelPseudopotentialModel::TermsFromFile(
        const io::fileFormat_t format,      //!<    Format of file
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t============ RETRIEVING TERM LOOK-UP TABLES FROM FILE ============"<<std::endl;
        }
        std::stringstream fileNameQuadratic, filenameQuartic;
        fileNameQuadratic.str("");
        filenameQuartic.str("");
        fileNameQuadratic << m_params.m_outPath << "/quadratic_coefficient_table_LLL_L_" << m_params.m_nbrOrbitals << ".dat";
        filenameQuartic << m_params.m_outPath << "/quartic_coefficient_table_LLL_L_" << m_params.m_nbrOrbitals << ".dat";
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            m_quadraticTables.FromFile(fileNameQuadratic.str(), format, mpi);
            m_quarticTables.FromFile(filenameQuartic.str(), format, mpi);
        }
        else
        {
            m_quadraticHashTables.FromFile(fileNameQuadratic.str(), format, mpi);
            m_quarticHashTables.FromFile(filenameQuartic.str(), format, mpi);
        }
        // 2LL terms
        fileNameQuadratic.str("");
        filenameQuartic.str("");
        fileNameQuadratic << m_params.m_outPath << "/quadratic_coefficient_table_2LL_L_" << m_params.m_nbrOrbitals << ".dat";
        filenameQuartic << m_params.m_outPath << "/quartic_coefficient_table_2LL_L_" << m_params.m_nbrOrbitals << ".dat";
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            m_quadraticTables2LL.FromFile(fileNameQuadratic.str(), format, mpi);
            m_quarticTables2LL.FromFile(filenameQuartic.str(), format, mpi);
        }
        else
        {
            m_quadraticHashTables2LL.FromFile(fileNameQuadratic.str(), format, mpi);
            m_quarticHashTables2LL.FromFile(filenameQuartic.str(), format, mpi);
        }
        MPI_Barrier(mpi.m_comm);
        m_params.m_termTablesBuilt = true;
        m_params.m_tableFormat = m_params.m_setTableFormat;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<" - DONE"<<std::endl;
        }
        return;
    }
    
}   //  End namespace diagonalization 
