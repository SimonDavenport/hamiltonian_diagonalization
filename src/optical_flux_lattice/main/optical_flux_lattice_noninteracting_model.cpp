////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This code performs some analysis of the non-interacting model
//!     constructed in the lowest lying band of an optical flux lattice. 
//!     See e.g. PRL 109, 265301 (2012).
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
#include "../../utilities/wrappers/program_options_wrapper.hpp"
#include "../noninteracting_ofl_model/noninteracting_ofl_model_grid.hpp"
#include "../../program_options/general_options.hpp"
#include "../program_options/sql_options.hpp"
///////     GLOBAL DATA STRUCTURES      ////////////////////////////////////////
utilities::Cout utilities::cout;
utilities::MpiWrapper mpi(utilities::cout);
///////		FUNCTION FORWARD DECLARATIONS		    ////////////////////////////
boost::program_options::variables_map ParseCommandLine(int argc, char *argv[],utilities::MpiWrapper& mpi);
///////		START OF MAIN FUNCTION      ////////////////////////////////////////
int main(int argc, char *argv[])
{
	mpi.Init(argc, argv);
    boost::program_options::variables_map optionList;
    optionList = ParseCommandLine(argc, argv, mpi);
    //////      GENERATE SINGLE PARTICLE MODEL DATA     ////////////////////////
    diagonalization::iSize_t spatialWaveFunctionsFlag = 0;
    bool calculateBandstructureFlag = false;
    bool calculateBandWidthFlag = false;
    bool calculateMagnetizationFlag = false;
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    GetOption(&optionList, spatialWaveFunctionsFlag, "calculate-spatial-wavefunctions", _AT_, mpi);
	    GetOption(&optionList, calculateBandstructureFlag, "calculate-bandstructure", _AT_, mpi);
	    GetOption(&optionList, calculateBandWidthFlag, "calculate-band-width", _AT_, mpi);
	    GetOption(&optionList, calculateMagnetizationFlag, "calculate-magnetization", _AT_, mpi);
	    bool plotBandStructureFlag;
	    bool plotBandWidthFlag;
	    bool plotMagnetizationFlag;
	    GetOption(&optionList, plotBandStructureFlag, "plot-bandstructure", _AT_, mpi);
	    GetOption(&optionList, plotBandWidthFlag, "plot-band-width", _AT_, mpi);
	    GetOption(&optionList, plotMagnetizationFlag, "plot-magnetization", _AT_, mpi);
	    calculateBandstructureFlag = calculateBandstructureFlag || plotBandStructureFlag;
	    calculateBandWidthFlag = calculateBandWidthFlag || plotBandWidthFlag;
	    calculateMagnetizationFlag = calculateMagnetizationFlag || plotMagnetizationFlag;
    }
    mpi.Sync(&spatialWaveFunctionsFlag, 1, 0);
    mpi.Sync(&calculateBandstructureFlag, 1, 0);
    mpi.Sync(&calculateBandWidthFlag, 1, 0);
    mpi.Sync(&calculateMagnetizationFlag, 1, 0);
    if(spatialWaveFunctionsFlag)
    {
        diagonalization::NonInteractingOflModelGrid modelGrid;
        modelGrid.CalculateSpatialWaveFunctions(&optionList, mpi);
    }
    if(calculateBandstructureFlag)
    {
        diagonalization::NonInteractingOflModelGrid modelGrid;
        modelGrid.PlotBandstructure(&optionList, mpi);
    }
    if(calculateBandWidthFlag)
    {
        diagonalization::NonInteractingOflModelGrid modelGrid;
        modelGrid.PlotBandWidth(&optionList, mpi);
    }
    if(calculateMagnetizationFlag)
    {
        diagonalization::NonInteractingOflModelGrid modelGrid;
        modelGrid.PlotMagnetization(&optionList, mpi);
    }
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
	    po::options_description noninteractingOflModelOpt(diagonalization::myOptions::GetCommonNonInteractingOflModelOptions());
	    diagonalization::myOptions::AddNonInteractingOflModelGridOptions(noninteractingOflModelOpt);
	    po::options_description allOpt("\n\tThis program generates the Hamiltonian describing a non-interacting optical flux lattice model.\n\tSee e.g. PRL109,265301 (2013)\n\n\tThe program input options are as follows");
	    allOpt.add(diagonalization::myOptions::GetGeneralOptions()).add(noninteractingOflModelOpt).add(diagonalization::myOptions::GetNonInteractingOflModelPlotOptions()).add(diagonalization::myOptions::GetCommonSqlOptions());
	    try
        {
            po::store(po::command_line_parser(argc,argv).options(allOpt).run(),vm);
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
	    GetOption(&vm, verbosity, "verbose", _AT_, mpi);
	    utilities::cout.SetVerbosity(verbosity);
    }
    utilities::cout.MpiSync(0, mpi.m_comm);
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        utilities::cout.MainOutput()<<"\n\tRun with -h option to see program options"<<std::endl;
    }
    return vm;
} 
