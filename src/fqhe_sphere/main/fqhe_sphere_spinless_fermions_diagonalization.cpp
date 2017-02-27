////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This code performs exact diagonalization of a FQHE pseudopotential
//!     Hamiltonian for spinless fermions in the sphere geometry
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
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../pseudopotential_hamiltonian/pseudopotential_hamiltonian.hpp"
#include "../pseudopotential_hamiltonian/two_level_pseudopotential_hamiltonian.hpp"
#include "../../program_options/general_options.hpp"
#include "../../program_options/arpack_options.hpp"
///////     GLOBAL DATA STRUCTURES      ////////////////////////////////////////
utilities::Cout utilities::cout;
utilities::MpiWrapper mpi(utilities::cout);
///////		FUNCTION FORWARD DECLARATIONS		    ////////////////////////////
boost::program_options::variables_map ParseCommandLine(int argc, char *argv[], 
                                                       utilities::MpiWrapper& mpi);
///////		START OF MAIN FUNCTION      ////////////////////////////////////////
int main(int argc, char *argv[])
{
	mpi.Init(argc, argv);
    boost::program_options::variables_map optionList;
    optionList = ParseCommandLine(argc, argv, mpi);
    //  Import and synchronize top level command line options over all nodes
    bool diagonalizeFlag;
    bool eigenvaluesFlag;
    bool eigenvectorsFlag;
    bool hamiltonianWriteFlag;
    diagonalization::iSize_t nbrLevels;
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    diagonalizeFlag  = optionList["diagonalize"].as<bool>();
	    eigenvaluesFlag  = optionList["eigenvalues-file"].as<bool>();
	    eigenvectorsFlag = optionList["eigenvectors-file"].as<bool>();
	    hamiltonianWriteFlag = optionList["hamiltonian-to-file"].as<bool>();
	    nbrLevels = optionList["nbr-levels"].as<diagonalization::iSize_t>();
	}
    //  MPI sync the flags from node 0
    mpi.Sync(&diagonalizeFlag, 1, 0);
    mpi.Sync(&eigenvaluesFlag, 1, 0);
    mpi.Sync(&eigenvectorsFlag, 1, 0);
    mpi.Sync(&hamiltonianWriteFlag, 1, 0);
    mpi.Sync(&nbrLevels, 1, 0);
    //////      BUILD AND DIAGONALIZE HAMILTONIAN       ////////////////////////
    if(diagonalizeFlag)
    {
        //  Determine angular momentum sectors to be diagonalized
        std::vector<diagonalization::iSize_t> sectorList;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {   
            if(optionList.count("lz-sectors"))
            {
                sectorList = optionList["lz-sectors"].as<std::vector<diagonalization::iSize_t> >();
            }
            else if(optionList["block-diagonalize"].as<bool>())
            {
                diagonalization::iSize_t startLz = 0;
                bool nbrParticlesOdd = optionList["nbr-particles"].as<diagonalization::iSize_t>() & 1;
                bool nbrOrbitalsOdd  = optionList["nbr-orbitals"].as<diagonalization::iSize_t>() & 1;
                if(!nbrOrbitalsOdd && nbrParticlesOdd)
                {
                    startLz = 1;
                }
                diagonalization::iSize_t nbrOrbitals =  optionList["nbr-orbitals"].as<diagonalization::iSize_t>();
                diagonalization::iSize_t nbrParticles = optionList["nbr-particles"].as<diagonalization::iSize_t>();
                diagonalization::iSize_t nbrLevels = optionList["nbr-levels"].as<diagonalization::iSize_t>();
                diagonalization::iSize_t max=0;
                if(nbrLevels==1)
                {
                    max = nbrParticles*(nbrOrbitals-nbrParticles);
                }
                else if(nbrLevels==2)
                {
                    max = (nbrOrbitals-2)/2+1;
                    for(unsigned int i=2; i<=nbrParticles; ++i)
                    {
                        max += ((nbrOrbitals-2)/2+1-2*(i/2));
                    }
                }
                for(diagonalization::iSize_t lz = startLz; lz<=max; lz+=1)
                {
                    sectorList.push_back(lz);
                }
            }
        }
        mpi.ExitFlagTest();
        mpi.Sync(&sectorList, 0);
        if(nbrLevels==1)
        {
            diagonalization::SpherePseudopotentialHamiltonian hamiltonian(&optionList, mpi);
            //  Generate look-up tables
            hamiltonian.BuildLookUpTables(mpi);
            if(0 == sectorList.size())
            {
                //  Construct the full Hamiltonian
                hamiltonian.BuildFockBasis(mpi);
                hamiltonian.BuildHamiltonian(mpi);
                if(hamiltonianWriteFlag)
                {
                    hamiltonian.HamiltonianToFile(mpi);
                }
                hamiltonian.Diagonalize(mpi);
                hamiltonian.EigensystemToFile(eigenvaluesFlag, eigenvectorsFlag, mpi);
            }
            else
            {
                //  Construct and diagonalize the Hamiltonian for each specified sector
                for(auto& sector : sectorList)
                {
                    //  Build and diagonalize selected sectors
                    hamiltonian.SetSector(sector);
                    hamiltonian.BuildFockBasis(mpi);
                    hamiltonian.BuildHamiltonian(mpi);
                    if(hamiltonianWriteFlag)
                    {
                        hamiltonian.HamiltonianToFile(mpi);
                    }
                    hamiltonian.Diagonalize(mpi);
                    hamiltonian.EigensystemToFile(eigenvaluesFlag, eigenvectorsFlag, mpi);
                    hamiltonian.ClearHamiltonian();
                }
            }
        }
        else if(nbrLevels==2)
        {
            diagonalization::SphereTwoLevelPseudopotentialHamiltonian hamiltonian(&optionList, mpi);
            //  Generate look-up tables
            hamiltonian.BuildLookUpTables(mpi);
            if(0 == sectorList.size())
            {
                //  Construct the full Hamiltonian
                hamiltonian.BuildHamiltonian(mpi);
                if(hamiltonianWriteFlag)
                {
                    hamiltonian.HamiltonianToFile(mpi);
                }
                hamiltonian.Diagonalize(mpi);
                hamiltonian.EigensystemToFile(eigenvaluesFlag, eigenvectorsFlag, mpi);
            }
            else
            {
                //  Construct and diagonalize the Hamiltonian for each specified sector
                for(auto& sector : sectorList)
                {
                    //  Build and diagonalize selected sectors
                    hamiltonian.SetSector(sector);
                    hamiltonian.BuildHamiltonian(mpi);
                    if(hamiltonianWriteFlag)
                    {
                        hamiltonian.HamiltonianToFile(mpi);
                    }
                    hamiltonian.Diagonalize(mpi);
                    hamiltonian.EigensystemToFile(eigenvaluesFlag, eigenvectorsFlag, mpi);
                    hamiltonian.ClearHamiltonian();
                }
            }
        }
    }  
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief A function to parse the command line arguments
//!
//! \return An instance of the boost program options variables map
//! containing the parsed command line arguments
////////////////////////////////////////////////////////////////////////////////
boost::program_options::variables_map ParseCommandLine(
    const int argc,                 //!<	Number of characters to parse
	char *argv[],	                //!<	Character array to parse
    utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
{
    namespace po = boost::program_options;
    po::variables_map vm;
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    //  Main program description
	    po::options_description allOpt("\n\tThis program generates a FQHE Haldane pseudopotential Hamiltonian for spinless fermions in the sphere geometry.\n\n\tThe program input options are as follows");
	    //	Declare option groups included
	    allOpt.add(diagonalization::myOptions::GetGeneralOptions()).add(diagonalization::myOptions::GetPseudopotentialHamiltonianOptions()).add(diagonalization::myOptions::GetArpackOptions());
	    try
        {
            po::store(po::command_line_parser(argc, argv).options(allOpt).run(), vm);
            if(vm.count("help"))
            {
	            utilities::cout.MainOutput()<<allOpt<<std::endl;
	            mpi.m_exitFlag = true;
            }
            po::notify(vm);
        }
        catch(po::error& e)
        {
            utilities::cout.MainOutput()<<allOpt<<std::endl;
            std::cerr<<utilities::cout.HyphenLine()<<std::endl;
            std::cerr<<std::endl<<"\tERROR:\t"<<e.what()<<std::endl;
            std::cerr<<std::endl<<utilities::cout.HyphenLine()<<std::endl;
            mpi.m_exitFlag = true;
        }
        //  Extract the options specifying the length of positional options and then re-parse with 
        //  the positional options included
    }
    mpi.ExitFlagTest();
    //  Set global verbosity level
	if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    utilities::cout.SetVerbosity(vm["verbose"].as<int>());
    }
    //  MPI sync verbosity level
    utilities::cout.MpiSync(0, mpi.m_comm);
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        utilities::cout.MainOutput()<<"\n\tRun with -h option to see program options"<<std::endl;
    }
    return vm;
}
