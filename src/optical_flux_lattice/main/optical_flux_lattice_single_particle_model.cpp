////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 31/10/2014
//!                                                                             
//!	 \file
//!     This code performs some analysis of the single particle model
//!     constructed in the lowest lying band of an optical flux lattice. 
//!     See e.g. PRL 109, 265301 (2012)
//!                                                                             
//!                    Copyright (C) 2014 Simon C Davenport
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

//	MPI functions and variables for parallelization
#include "../../utilities/wrappers/mpi_wrapper.hpp"

//  Functions to manipulate std::cout output
#include "../../utilities/general/cout_tools.hpp"

//  A class to contain an array of single particle models (at different points
//  on a k-space grid)
#include "../single_particle_hamiltonian/single_particle_hamiltonian_array.hpp"

//  Include program options declarations
#include "../../program_options/general_options.hpp"
#include "../program_options/sql_options.hpp"

///////     GLOBAL DATA STRUCTURES      ////////////////////////////////////////

//  Declare an instance of the global utilities::cout struct for cout_tools
utilities::Cout utilities::cout;

//  Declare an instance of the global mpi wrapper class
utilities::MpiWrapper mpi(utilities::cout);

///////		FUNCTION FORWARD DECLARATIONS		    ////////////////////////////

//  Declare a function to define the command line arguments
boost::program_options::variables_map ParseCommandLine(int argc, char *argv[],utilities::MpiWrapper& mpi);

///////		START OF MAIN FUNCTION      ////////////////////////////////////////

int main(int argc, char *argv[])
{
    //////  	START PARALLEL PROCESSES      //////////////////////////////////

	mpi.Init(argc,argv);
	
	//  Declare a command line variables map on each node
    boost::program_options::variables_map optionList;

    //////      PARSE COMMAND LINE ARGUMENTS    //////

    optionList = ParseCommandLine(argc,argv,mpi);
    
    //////      GENERATE SINGLE PARTICLE MODEL DATA     ////////////////////////
    
    diagonalization::iSize_t spatialWaveFunctionsFlag = 0;
    bool calculateBandstructureFlag = false;
    bool calculateBandWidthFlag     = false;
    bool calculateMagnetizationFlag = false;
    
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    spatialWaveFunctionsFlag   = optionList["calculate-spatial-wavefunctions"].as<diagonalization::iSize_t>();
        calculateBandstructureFlag = optionList["calculate-bandstructure"].as<bool>() || optionList["plot-bandstructure"].as<bool>();
        calculateBandWidthFlag     = optionList["calculate-band-width"].as<bool>() || optionList["plot-band-width"].as<bool>();
        calculateMagnetizationFlag = optionList["calculate-magnetization"].as<bool>() || optionList["plot-magnetization"].as<bool>();
    }
    
    mpi.Sync(&spatialWaveFunctionsFlag,1,0);
    mpi.Sync(&calculateBandstructureFlag,1,0);
    mpi.Sync(&calculateBandWidthFlag,1,0);
    mpi.Sync(&calculateMagnetizationFlag,1,0);
    
    if(spatialWaveFunctionsFlag)
    {
        diagonalization::SingleParticleHamiltonianArray modelArray;
    
        modelArray.CalculateSpatialWaveFunctions(&optionList,mpi);
    }
    
    if(calculateBandstructureFlag)
    {
        diagonalization::SingleParticleHamiltonianArray modelArray;
        
        modelArray.PlotBandstructure(&optionList,mpi);
    }
    
    if(calculateBandWidthFlag)
    {
        diagonalization::SingleParticleHamiltonianArray modelArray;
    
        modelArray.PlotBandWidth(&optionList,mpi);
    }
    
    if(calculateMagnetizationFlag)
    {
        diagonalization::SingleParticleHamiltonianArray modelArray;
    
        modelArray.PlotMagnetization(&optionList,mpi);
    }
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief A function to parse the command line arguments
//!
//! \return An instance of the boost program options variables map
//! containing the parsed command line arguments
//!
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
	    //  Construct single particle model options list
	    po::options_description sinlgeParticleOpt(diagonalization::myOptions::GetCommonSingleParticleOptions());
	    
	    diagonalization::myOptions::AddSingleParticleArrayOptions(sinlgeParticleOpt);
	
	    //  Main program description
	    po::options_description allOpt("\n\tThis program generates the single particle Hamiltonian describing an optical flux lattice model.\n\tSee e.g. PRL109,265301 (2013)\n\n\tThe program input options are as follows");

	    //	Declare option groups included
	    allOpt.add(diagonalization::myOptions::GetGeneralOptions()).add(sinlgeParticleOpt).add(diagonalization::myOptions::GetSingleParticlePlotOptions()).add(diagonalization::myOptions::GetCommonSqlOptions());

	    try
        {
            //	Map the command line arguments
            po::store(po::command_line_parser(argc,argv).options(allOpt).run(),vm);
            
            //	Respond to help request
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
    
    //  Set global verbosity level
	
	if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    utilities::cout.SetVerbosity(vm["verbose"].as<int>());
    }
    
    //  MPI sync verbosity level
    
    utilities::cout.MpiSync(0,mpi.m_comm);

    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        utilities::cout.MainOutput()<<"\n\tRun with -h option to see program options"<<std::endl;
    }

    return vm;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//    
