////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This code performs exact diagonalization of an interacting model
//!     constructed in the lowest lying band of an optical flux lattice.
//!     See e.g. PRL 109, 265301 (2012)
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
#include <vector>
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../interacting_ofl_model/interacting_ofl_model.hpp"
#include "../../program_options/general_options.hpp"
#include "../program_options/sql_options.hpp"
#include "../observables_manager/observables_manager.hpp"
#include "../../utilities/wrappers/program_options_wrapper.hpp"
///////     GLOBAL DATA STRUCTURES      ////////////////////////////////////////
utilities::Cout utilities::cout;
utilities::MpiWrapper mpi(utilities::cout);
///////		FUNCTION FORWARD DECLARATIONS		    ////////////////////////////
boost::program_options::variables_map ParseCommandLine(int argc, char *argv[], utilities::MpiWrapper& mpi);
std::vector<std::complex<diagonalization::iSize_t> > GenerateSectorList(boost::program_options::variables_map* optionList,
                                                                        utilities::MpiWrapper& mpi);
///////		START OF MAIN FUNCTION      ////////////////////////////////////////
int main(int argc, char *argv[])
{
	mpi.Init(argc, argv);
    boost::program_options::variables_map optionList;
    optionList = ParseCommandLine(argc, argv, mpi);
    //  Import and synchronize top-level command line options over all nodes
    bool diagonalizeFlag;
    bool eigenvaluesFlag;
    bool eigenvectorsFlag;
    bool storeTermsFlag;
    bool retrieveTermsFlag;
    int formatCode;
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    GetOption(&optionList, diagonalizeFlag, "diagonalize", _LINE_, mpi);
	    GetOption(&optionList, eigenvaluesFlag, "eigenvalues-file", _LINE_, mpi);
	    GetOption(&optionList, eigenvectorsFlag, "eigenvectors-file", _LINE_, mpi);
	    GetOption(&optionList, storeTermsFlag, "store-terms", _LINE_, mpi);
	    GetOption(&optionList, retrieveTermsFlag, "retrieve-terms", _LINE_, mpi);
	    GetOption(&optionList, formatCode, "file-format", _LINE_, mpi);
	}
	mpi.ExitFlagTest();
    mpi.Sync(&diagonalizeFlag, 1, 0);
    mpi.Sync(&eigenvaluesFlag, 1, 0);
    mpi.Sync(&eigenvectorsFlag, 1, 0);
    mpi.Sync(&storeTermsFlag, 1, 0);
    mpi.Sync(&retrieveTermsFlag, 1, 0);
    mpi.Sync(&formatCode, 1, 0);
    io::fileFormat_t fileFormat = diagonalization::myOptions::GetFileFormat(formatCode, mpi);
    //  Make a list of observables to calculate
    diagonalization::ObservablesManager observables;
    observables.AddObservables(&optionList, mpi);
    if(observables.GetNbrSetObservables()>0)
    {
        observables.Print(mpi);
    }
    //////      BUILD AND DIAGONALIZE HAMILTONIAN       ////////////////////////
    if(diagonalizeFlag)
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.MainOutput() << "\n\t============ DIAGONALIZATION ============" << std::endl;
        }
        std::vector<std::complex<diagonalization::iSize_t> > sectorList = GenerateSectorList(&optionList, mpi);
        diagonalization::InteractingOflModel model(&optionList, mpi);
        if(retrieveTermsFlag)
        {
            model.TermsFromFile(fileFormat, mpi);
        }
        else
        {
            model.BuildTermTables(&optionList, mpi);
        }
        if(0 == sectorList.size())  //  Diagonalize all sectors
        {
            model.BuildFockBasis(mpi);
            model.BuildHamiltonian(mpi);
            model.Diagonalize(mpi);
            model.EigensystemToFile(eigenvaluesFlag, eigenvectorsFlag, fileFormat, mpi);
            observables.CalculateAllObservables(&model, mpi);
        }
        else    //  Diagonalize selected sectors
        {
            for(std::vector<std::complex<diagonalization::iSize_t> >::const_iterator it = sectorList.begin(); 
            it != sectorList.end(); ++it)
            {
                model.SetSector(it->real(),it->imag());
                model.BuildFockBasis(mpi);
                model.BuildHamiltonian(mpi);
                model.Diagonalize(mpi);
                model.EigensystemToFile(eigenvaluesFlag, eigenvectorsFlag, fileFormat, mpi);
                observables.CalculateAllObservables(&model, mpi);
                model.ClearHamiltonian();
            }
        }
        model.UpdateSqlStatus(mpi);
        observables.UpdateSqlFlags(&optionList, mpi);
    }
    //////      COMPUTE OBSERVABLES USING STORED EIGENSTATES       /////////////
    else if(observables.GetNbrSetObservables()>0)   // Use existing results for analysis
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        utilities::cout.MainOutput() << "\n\t============ COMPUTE OBSERVABLES ============" << std::endl;
            utilities::cout.MainOutput()<<"\n\tLooking for existing eigensystem data to analyse."<<std::endl;
        }
        std::vector<std::complex<diagonalization::iSize_t> > sectorList = GenerateSectorList(&optionList, mpi);
        diagonalization::InteractingOflModel model(&optionList, mpi);
        bool calculateSusceptibilityFlag = false;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        GetOption(&optionList, calculateSusceptibilityFlag, "calculate-susceptibility", _LINE_, mpi);
	    }
        mpi.Sync(&calculateSusceptibilityFlag, 1, 0);
        if(calculateSusceptibilityFlag)
        {
            if(retrieveTermsFlag)
            {
                model.TermsFromFile(fileFormat, mpi);
            }
            else
            {
                model.BuildTermTables(&optionList, mpi);
            }
        }
        eigenvaluesFlag = true;
        eigenvectorsFlag = true;
        if(0 == sectorList.size())
        {
            model.EigensystemFromFile(eigenvaluesFlag, eigenvectorsFlag, fileFormat, mpi);
            model.BuildFockBasis(mpi);
            observables.CalculateAllObservables(&model, mpi);
        }
        else
        {
            for(std::vector<std::complex<diagonalization::iSize_t> >::const_iterator it = sectorList.begin(); 
            it != sectorList.end(); ++it)
            {
                model.SetSector(it->real(),it->imag());
                model.EigensystemFromFile(eigenvaluesFlag, eigenvectorsFlag, fileFormat, mpi);
                model.BuildFockBasis(mpi);
                observables.CalculateAllObservables(&model, mpi);
            }
        }
        model.UpdateSqlStatus(mpi);
        observables.UpdateSqlFlags(&optionList, mpi);
    }
    //////      GENERATE AND STORE TERM TABLES       ///////////////////////////
    if(storeTermsFlag)
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.MainOutput() << "\n\t============ COMPUTE TERMS FOR STORAGE ============" << std::endl;
        }
        diagonalization::InteractingOflModel model(&optionList, mpi);
        model.BuildTermTables(&optionList, mpi);
        model.TermsToFile(fileFormat, mpi);
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////
//! \brief A function to parse the command line arguments
//!
//! \return An instance of the boost program options variables map
//! containing the parsed command line arguments
//////////////////////////////////////////////////////////////////////////////////
boost::program_options::variables_map ParseCommandLine(
    const int argc,             //!<	Number of characters to parse
	char *argv[],	            //!<	Character array to parse
	utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
{
    namespace po = boost::program_options;
    po::variables_map vm;
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        po::options_description modelOpt(diagonalization::myOptions::GetCommonInteractingModelOptions());
        diagonalization::myOptions::AddTermOptions(modelOpt);
        diagonalization::myOptions::AddMoreInteractingModelOptions(modelOpt);
        po::options_description sqlOpt(diagonalization::myOptions::GetCommonSqlOptions());
        po::options_description allOpt("\n\tThis program generates the interacting Hamiltonian due to k-dependent interactions\n\tin the lowest lying band of an optical flux lattice model.\n\tSee e.g. PRL109,265301 (2013)\n\n\tThe program input options are as follows");
        allOpt.add(diagonalization::myOptions::GetGeneralOptions()).add(diagonalization::myOptions::GetCommonNonInteractingOflModelOptions()).add(modelOpt).add(diagonalization::myOptions::GetInteractingModelPlotOptions()).add(diagonalization::myOptions::GetObservablesOptions()).add(sqlOpt).add(diagonalization::myOptions::GetArpackOptions());
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
    }
    mpi.ExitFlagTest();
	if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    int verbosity;
	    GetOption(&vm, verbosity, "verbose", _LINE_, mpi);
	    utilities::cout.SetVerbosity(verbosity);
    }
    utilities::cout.MpiSync(0, mpi.m_comm);
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        utilities::cout.MainOutput()<<"\n\tRun with -h option to see program options"<<std::endl;
    }
    return vm;
}

//////////////////////////////////////////////////////////////////////////////////
//! \brief A function to generate a list of pairs kx_tot ky_tot  (or just ky_tot
//! in the Wannier basis) of linear momentum sectors to be diagonalized
//!
//! \return a list of pairs of kx_tot ky_tot  (or just ky_tot) linear momentum sectors
//////////////////////////////////////////////////////////////////////////////////
std::vector<std::complex<diagonalization::iSize_t> > GenerateSectorList(
    boost::program_options::variables_map* optionList,
                                //!<    Import command line option list to set sectors
    utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
{
    diagonalization::iSize_t* sectorList      = 0;
    diagonalization::iSize_t lengthSectorList = 0;
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    std::vector<diagonalization::iSize_t> tempList;
	    bool blockDiagonalize;
	    GetOption(optionList, blockDiagonalize, "block-diagonalize", _LINE_, mpi);
	    if(optionList->count("sectors"))
	    {
	        utilities::GetOption(optionList, tempList, "sectors", _LINE_, mpi);
        }
        else if(blockDiagonalize)
        {
            diagonalization::iSize_t kx;
            diagonalization::iSize_t ky;
            int basisCode;
            GetOption(optionList, kx, "kx", _LINE_, mpi);
            GetOption(optionList, ky, "ky", _LINE_, mpi);
            GetOption(optionList, basisCode, "basis", _LINE_, mpi);
            if(0 == basisCode)  // use kx,ky basis
            {
                for(diagonalization::iSize_t x=0;x<kx;++x)
                {
                    for(diagonalization::iSize_t y=0;y<ky;++y)
                    {
                        tempList.push_back(x);
                        tempList.push_back(y);
                    }
                }
            }
            else if(1 == basisCode) //  use ky,x basis
            {
                diagonalization::iSize_t x = 0;
                for(diagonalization::iSize_t y=0;y<ky;++y)
                {
                    tempList.push_back(x);
                    tempList.push_back(y);
                }
            }
        }
        if((tempList.size() & 1) != 0)
        {
            std::cerr<<"\n\tERROR WITH \"sectors\" ARGUMENT: odd number of k values specified."
                <<"\n\tExpecting pairs of kx_tot ky_tot [kx_tot ignored for Wannier basis case]."<<std::endl;
            mpi.m_exitFlag = true;
        }
        //  Convert to int array (this is required in order to use MPI functions)
        sectorList = new diagonalization::iSize_t[tempList.size()];
        for(diagonalization::iSize_t i=0; i<tempList.size(); ++i)
        {
            sectorList[i] = tempList[i];
        }
        lengthSectorList = tempList.size();
    }
    mpi.ExitFlagTest();
    mpi.Sync(&lengthSectorList, 1, 0);
    if(0 != mpi.m_id)	// FOR THE MASTER NODE
	{
	    sectorList = new diagonalization::iSize_t[lengthSectorList];
	}
    mpi.Sync(sectorList,lengthSectorList, 0);
    //  Finally convert to a std::vector<std::complex<iSize_t> > for return
    std::vector<std::complex<diagonalization::iSize_t> > returnArray;
    for(diagonalization::iSize_t i=0; i<lengthSectorList/2; ++i)
    {
        returnArray.push_back(std::complex<diagonalization::iSize_t>(sectorList[2*i], sectorList[2*i+1]));
    }
    delete[] sectorList;
    return returnArray;
}
