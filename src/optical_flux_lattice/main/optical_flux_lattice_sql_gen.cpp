////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 06/11/2014
//!                                                                             
//!	 \file
//!     This code generates an SQL database that is designed to keep track of
//!     a large set of diagonalizations performed by the 
//!     optical_flux_lattice_diag program                         
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

//  Sqlite wrapper functions for SQL data management
#include "../../utilities/wrappers/sqlite_wrapper.hpp"

//  Functions to manipulate std::cout output
#include "../../utilities/general/cout_tools.hpp"

//  For printing out a load bar
#include "../../utilities/general/load_bar.hpp"

//  Include program options declarations
#include "../../program_options/general_options.hpp"
#include "../program_options/sql_options.hpp"
#include "../program_options/optical_flux_lattice_hamiltonian_options.hpp"

///////     GLOBAL DATA STRUCTURES      ////////////////////////////////////////

//  Declare an instance of the global utilities::cout struct for cout_tools
utilities::Cout utilities::cout;

///////		FUNCTION FORWARD DECLARATIONS		    ////////////////////////////

//  Declare a function to define the command line arguments
boost::program_options::variables_map ParseCommandLine(int argc, char *argv[]);

///////		START OF MAIN FUNCTION      ////////////////////////////////////////

int main(int argc, char *argv[])
{
    //  Declare a command line variables map on each node
    boost::program_options::variables_map optionList;

    //////      PARSE COMMAND LINE ARGUMENTS    //////

    optionList = ParseCommandLine(argc,argv);

    if(optionList["build-sql-table"].as<bool>())
    {
        //////      GENERATE AN SQL TABLE FOR A CUT THROUGH PARAMETER SPACE     ////////
     
        utilities::cout.MainOutput()<<"\n\t============ BUILDING SQL TABLE ==========\n"<<std::endl;

        //////      SET PARAMETER CUT VALUES HERE       ////////////

        //  Grid of lattice depth vs interaction strength
        
        const diagonalization::iSize_t xGrid  = optionList["sql-v0-nbr"].as<diagonalization::iSize_t>();
        const diagonalization::iSize_t yGrid  = optionList["sql-g-nbr"].as<diagonalization::iSize_t>();
        const double v0Min           = optionList["sql-v0-min"].as<double>();
        const double v0Step          = optionList["sql-v0-step"].as<double>();
        const double interactionMin  = optionList["sql-g-min"].as<double>();
        const double interactionStep = optionList["sql-g-step"].as<double>();
        const double offsetX         = optionList["kx-shift"].as<double>();
        const double offsetY         = optionList["ky-shift"].as<double>();
        
        //  Other parameters are fixed
        
        const double defaultTheta   = 0.3;
        const double defaultEpsilon = 0.4;
        const double defaultKappa   = 1;
        const double defaultMass    = 1;

        utilities::cout.SecondaryOutput()<<"\tV0 from "<<v0Min<<" to "<<(v0Min+v0Step*(xGrid-1))<<" in steps of "<<v0Step<<std::endl;
        utilities::cout.SecondaryOutput()<<"\ttilde g/2Pi from "<<interactionMin<<" to "<<(interactionMin+interactionStep*(yGrid-1))<<" in steps of "<<interactionStep<<std::endl;
        utilities::cout.SecondaryOutput()<<"\ttheta = "<<defaultTheta<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tepsilon = "<<defaultEpsilon<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tkappa = "<<defaultKappa<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tmass = "<<defaultMass<<std::endl;

        utilities::cout.SecondaryOutput()<<"\n\tGENERATING SQL TABLE ";

        //////      GENERATE SQL FILE       ////////////
        
        std::stringstream fileName;
        fileName.str("");

        //  Generate an sql database with the following name
        fileName<<optionList["out-path"].as<std::string>()<<optionList["sql-name"].as<std::string>();

        utilities::cout.SecondaryOutput()<<fileName.str()<<std::endl<<std::endl;

        //  the _CREATE_NEW_ flag either creates a new sql file, or appends to existing file
        
        utilities::Sqlite sql(fileName.str(),_CREATE_NEW_);
        
        sql.AddTimeUpdatedTrigger("TimeEntered");

        //////      DEFINE TABLE DATA FIELDS        //////
        
        utilities::SqliteRow sqlRow;

        sqlRow.AddField("Theta",defaultTheta);       //  Theta 
        sqlRow.AddField("V0",_REAL_);                //  Lattice depth
        sqlRow.AddField("Epsilon",defaultEpsilon);   //  Epsilon 
        sqlRow.AddField("Kappa",defaultKappa);       //  Kappa
        sqlRow.AddField("Mass",defaultMass);         //  Mass
        sqlRow.AddField("Interaction",_REAL_);       //  Interaction strength
        sqlRow.AddField("OutFileName",_TEXT_);       //  Output file name
        sqlRow.AddField("OffsetX",offsetX);          //  kx offset
        sqlRow.AddField("OffsetY",offsetY);          //  ky offset
        sqlRow.AddField("JobId",_INT_);              //  Leave a blank column to record the associated cluster job id
        
        int completedFlag = 0;
        
        if(optionList["sql-completed"].as<bool>())
        {
            completedFlag = 1;
        }
        
        //  Flag to say whether the task is completed
        
        sqlRow.AddField("Completed",completedFlag);
     
        //  Flags for additional data analysis
        
        sqlRow.AddField("GotOccupations",completedFlag);
        sqlRow.AddField("GotSusceptibility",completedFlag);
        //sqlRow.AddField("GotOverlaps",completedFlag);
        sqlRow.AddField("GotMostProbable",completedFlag);
        sqlRow.AddField("GotDensityDensity",completedFlag);
        sqlRow.AddField("GotParticipationRatio",completedFlag);
        //sqlRow.AddField("GotOES",completedFlag);
        //sqlRow.AddField("GotPES",completedFlag);
        
        //  Generate a table to contain these data
      
        std::string tableName = optionList["sql-table-name"].as<std::string>();
        
        sql.CreateTable(tableName,&sqlRow);

        //////      POPULATE THE TABLE WITH A CUT THROUGH PARAMETER SPACE   ///////

        double currV0 = v0Min;
        diagonalization::iSize_t startCounter = sql.GetMaxId(tableName)+1;
        diagonalization::iSize_t idCounter    = startCounter;
        
        utilities::LoadBar progress;
        
        progress.Initialize(xGrid*yGrid);
        
        //  Create a zig-zag path through parameter space to keep changes in parameter 
        //  values minimal (thus allowing an accurate analysis of eigenvector overlap
        //  matrices)
        
        //  A list of all data to be entered into the table

        std::vector<utilities::SqliteRow> data;
        
        data.reserve(xGrid*yGrid);

        for(diagonalization::iSize_t x=0;x<xGrid;++x,currV0+=v0Step)
        {
            sqlRow.UpdateValue("V0",currV0);

            if(x&1) //  Snake one way if odd
            {
                double currInteraction = (interactionMin+interactionStep*(yGrid-1));
            
                for(diagonalization::iSize_t y=0;y<yGrid;y++,currInteraction-=interactionStep,++idCounter)
                {
                    std::stringstream ss;
                
                    ss.str("");

                    ss<<"optical_flux_model_n_"<<optionList["nbr"].as<diagonalization::iSize_t>()
                    <<"_kx_"<<optionList["kx"].as<diagonalization::iSize_t>()<<"_ky_"<<optionList["ky"].as<diagonalization::iSize_t>()
                    <<"_id_"<<idCounter;
                
                    if(fabs(currInteraction)<10e-15)    currInteraction = 0.0;
                
                    sqlRow.UpdateValue("Interaction",currInteraction);
                    
                    sqlRow.UpdateValue("OutFileName",_TEXT_,ss.str());
                    
                    data.push_back(sqlRow);

                    progress.Display(idCounter-startCounter+1);
                }
            }
            else    //  Snake the other way if even
            {
                double currInteraction = interactionMin;
            
                for(diagonalization::iSize_t y=0;y<yGrid;++y,currInteraction+=interactionStep,++idCounter)
                {
                    std::stringstream ss;
                
                    ss.str("");

                    ss<<"optical_flux_model_n_"<<optionList["nbr"].as<diagonalization::iSize_t>()
                    <<"_kx_"<<optionList["kx"].as<diagonalization::iSize_t>()<<"_ky_"<<optionList["ky"].as<diagonalization::iSize_t>()
                    <<"_id_"<<idCounter;
                
                    if(fabs(currInteraction)<10e-15)    currInteraction = 0.0;
                
                    sqlRow.UpdateValue("Interaction",currInteraction);
                    
                    sqlRow.UpdateValue("OutFileName",_TEXT_,ss.str());
                    
                    data.push_back(sqlRow);
     
                    progress.Display(idCounter-startCounter+1);
                }
            }
        }
        
        //  Insert all data into the SQL table
        
        sql.InsertIntoTable(tableName,&data);

        sql.CopyTableToTextFile(tableName,fileName.str());
        
        utilities::cout.MainOutput()<<"\n\tDONE\n"<<std::endl;
    }
    else if(optionList["build-sql-table-offset"].as<bool>())
    {
        //////      GENERATE AN SQL TABLE FOR A CUT THROUGH OFFSET SPACE     ////////
     
        utilities::cout.MainOutput()<<"\n\t============ BUILDING SQL TABLE ==========\n"<<std::endl;
        
        //////      SET PARAMETER CUT VALUES HERE       ////////////

        const double defaultV0       = optionList["sql-v0-min"].as<double>();
        const double defaultG        = optionList["sql-g-min"].as<double>();
        const double offsetXmin      = optionList["sql-kx-shift-min"].as<double>();
        const diagonalization::iSize_t offsetXnbr = optionList["sql-kx-shift-nbr"].as<diagonalization::iSize_t>();
        const double offsetXstep     = optionList["sql-kx-shift-step"].as<double>();
        const double offsetYmin      = optionList["sql-ky-shift-min"].as<double>();
        const diagonalization::iSize_t offsetYnbr = optionList["sql-ky-shift-nbr"].as<diagonalization::iSize_t>();
        const double offsetYstep     = optionList["sql-ky-shift-step"].as<double>();
        
        //  Other parameters are fixed
        
        const double defaultTheta   = 0.3;
        const double defaultEpsilon = 0.4;
        const double defaultKappa   = 1;
        const double defaultMass    = 1;

        utilities::cout.SecondaryOutput()<<"\tkx offset from "<<offsetXmin<<" to "<<(offsetXmin+offsetXstep*(offsetXnbr-1))<<" in steps of "<<offsetXstep<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tky offset from "<<offsetYmin<<" to "<<(offsetYmin+offsetYstep*(offsetYnbr-1))<<" in steps of "<<offsetYstep<<std::endl;  
        utilities::cout.SecondaryOutput()<<"\tV0 = "<<defaultV0<<std::endl;
        utilities::cout.SecondaryOutput()<<"\ttilde g/2Pi = "<<defaultG<<std::endl;
        utilities::cout.SecondaryOutput()<<"\ttheta = "<<defaultTheta<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tepsilon = "<<defaultEpsilon<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tkappa = "<<defaultKappa<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tmass = "<<defaultMass<<std::endl;

        utilities::cout.SecondaryOutput()<<"\n\tGENERATING SQL TABLE ";
        
        //////      GENERATE SQL FILE       ////////////
        
        std::stringstream fileName;
        fileName.str("");

        //  Generate an sql database with the following name
        fileName<<optionList["out-path"].as<std::string>()<<optionList["sql-name"].as<std::string>();

        utilities::cout.SecondaryOutput()<<fileName.str()<<std::endl<<std::endl;

        //  the _CREATE_NEW_ flag either creates a new sql file, or appends to existing file
        
        utilities::Sqlite sql(fileName.str(),_CREATE_NEW_);
        
        sql.AddTimeUpdatedTrigger("TimeEntered");

        //////      DEFINE TABLE DATA FIELDS        //////
        
        utilities::SqliteRow sqlRow;

        sqlRow.AddField("Theta",defaultTheta);       //  Theta 
        sqlRow.AddField("V0",defaultV0);             //  Lattice depth
        sqlRow.AddField("Epsilon",defaultEpsilon);   //  Epsilon 
        sqlRow.AddField("Kappa",defaultKappa);       //  Kappa
        sqlRow.AddField("Mass",defaultMass);         //  Mass
        sqlRow.AddField("Interaction",defaultG);     //  Interaction strength
        sqlRow.AddField("OutFileName",_TEXT_);       //  Output file name
        sqlRow.AddField("OffsetX",_REAL_);           //  kx offset
        sqlRow.AddField("OffsetY",_REAL_);           //  ky offset
        sqlRow.AddField("JobId",_INT_);              //  Leave a blank column to record the associated cluster job id
        
        int completedFlag = 0;
        
        if(optionList["sql-completed"].as<bool>())
        {
            completedFlag = 1;
        }
        
        //  Flag to say whether the task is completed
        
        sqlRow.AddField("Completed",completedFlag);
     
        //  Flags for additional data analysis
        
        sqlRow.AddField("GotOccupations",completedFlag);
        sqlRow.AddField("GotSusceptibility",completedFlag);
        //sqlRow.AddField("GotOverlaps",completedFlag);
        sqlRow.AddField("GotMostProbable",completedFlag);
        sqlRow.AddField("GotDensityDensity",completedFlag);
        sqlRow.AddField("GotParticipationRatio",completedFlag);
        //sqlRow.AddField("GotOES",completedFlag);
        //sqlRow.AddField("GotPES",completedFlag);
        
        //  Generate a table to contain these data
      
        std::string tableName = optionList["sql-table-name"].as<std::string>();
        
        sql.CreateTable(tableName,&sqlRow);

        //////      POPULATE THE TABLE WITH A CUT THROUGH OFFSET SPACE   ///////
        
        double currxOffsetX = offsetXmin;
        diagonalization::iSize_t startCounter = sql.GetMaxId(tableName)+1;
        diagonalization::iSize_t idCounter    = startCounter;
        
        utilities::LoadBar progress;
        
        progress.Initialize(offsetXnbr*offsetYnbr);

        //  A list of all data to be entered into the table

        std::vector<utilities::SqliteRow> data;
        
        data.reserve(offsetXnbr*offsetYnbr);

        for(diagonalization::iSize_t x=0;x<offsetXnbr;++x,currxOffsetX+=offsetXstep)
        {
            sqlRow.UpdateValue("OffsetX",currxOffsetX);
            
            double currxOffsetY = offsetYmin;
            
            for(diagonalization::iSize_t y=0;y<offsetYnbr;++y,currxOffsetY+=offsetYstep,++idCounter)
            {
                std::stringstream ss;
                
                ss.str("");

                ss<<"optical_flux_model_n_"<<optionList["nbr"].as<diagonalization::iSize_t>()
                <<"_kx_"<<optionList["kx"].as<diagonalization::iSize_t>()<<"_ky_"<<optionList["ky"].as<diagonalization::iSize_t>()
                <<"_id_"<<idCounter;

                sqlRow.UpdateValue("OffsetY",currxOffsetY);
                
                sqlRow.UpdateValue("OutFileName",_TEXT_,ss.str());
                
                data.push_back(sqlRow);

                progress.Display(idCounter-startCounter+1);
            }
        }
        
        //  Insert all data into the SQL table
        
        sql.InsertIntoTable(tableName,&data);

        sql.CopyTableToTextFile(tableName,fileName.str());
        
        utilities::cout.MainOutput()<<"\n\tDONE\n"<<std::endl;
    }
    else if(optionList["build-sql-table-single-particle"].as<bool>())
    {
        //////      GENERATE AN SQL TABLE FOR A CUT THROUGH PARAMETER SPACE     ////////
     
        utilities::cout.MainOutput()<<"\n\t============ BUILDING SQL TABLE ==========\n"<<std::endl;

        //////      SET PARAMETER CUT VALUES HERE       ////////////
        
        const double offsetX         = optionList["kx-shift"].as<double>();
        const double offsetY         = optionList["ky-shift"].as<double>();

        const diagonalization::iSize_t epsilonNbr = 10;
        const double epsilonMin          = 0.1;
        const double epsilonStep         = 0.1;
        
        const diagonalization::iSize_t thetaNbr  = 10;
        const double thetaMin           = 0.1;
        const double thetaStep          = 0.1;
        const double defaultV0          = 1.0;
        
        //  Other parameters are fixed

        const double defaultKappa   = 1;
        const double defaultMass    = 1;

        utilities::cout.SecondaryOutput()<<"\tEpsilon from "<<epsilonMin<<" to "<<(epsilonMin+epsilonStep*(epsilonNbr-1))<<" in steps of "<<epsilonStep<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tTheta from "<<thetaMin<<" to "<<(thetaMin+epsilonStep*(thetaNbr-1))<<" in steps of "<<thetaStep<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tkappa = "<<defaultKappa<<std::endl;
        utilities::cout.SecondaryOutput()<<"\tmass = "<<defaultMass<<std::endl;

        utilities::cout.SecondaryOutput()<<"\n\tGENERATING SQL TABLE ";
        
        //////      GENERATE SQL FILE       ////////////
        
        std::stringstream fileName;
        fileName.str("");

        //  Generate an sql database with the following name
        fileName<<optionList["out-path"].as<std::string>()<<optionList["sql-name"].as<std::string>();

        utilities::cout.SecondaryOutput()<<fileName.str()<<std::endl<<std::endl;

        //  the _CREATE_NEW_ flag either creates a new sql file, or appends to existing file
        
        utilities::Sqlite sql(fileName.str(),_CREATE_NEW_);
        
        sql.AddTimeUpdatedTrigger("TimeEntered");

        //////      DEFINE TABLE DATA FIELDS        //////
        
        utilities::SqliteRow sqlRow;

        sqlRow.AddField("Epsilon",_REAL_);           //  Epsilon
        sqlRow.AddField("Theta",_REAL_);             //  Theta 
        sqlRow.AddField("V0",defaultV0);             //  Lattice depth
        sqlRow.AddField("Kappa",defaultKappa);       //  Kappa
        sqlRow.AddField("Mass",defaultMass);         //  Mass
        sqlRow.AddField("OutFileName",_TEXT_);       //  Output file name
        sqlRow.AddField("OffsetX",offsetX);          //  kx offset
        sqlRow.AddField("OffsetY",offsetY);          //  ky offset
        sqlRow.AddField("JobId",_INT_);              //  Leave a blank column to record the associated cluster job id
        
        int completedFlag = 0;
        
        if(optionList["sql-completed"].as<bool>())
        {
            completedFlag = 1;
        }
        
        //  Flag to say whether the task is completed
        
        sqlRow.AddField("Completed",completedFlag);
        
        //  Generate a table to contain these data
      
        std::string tableName = optionList["sql-table-name"].as<std::string>();
        
        sql.CreateTable(tableName,&sqlRow);

        //////      POPULATE THE TABLE WITH A CUT THROUGH PARAMETER SPACE   ///////

        double currEpsilon = epsilonMin;
        diagonalization::iSize_t startCounter = sql.GetMaxId(tableName)+1;
        diagonalization::iSize_t idCounter    = startCounter;
        
        utilities::LoadBar progress;
        
        progress.Initialize(epsilonNbr*thetaNbr);
        
        //  A list of all data to be entered into the table

        std::vector<utilities::SqliteRow> data;
        
        data.reserve(epsilonNbr*thetaNbr);
        
        for(diagonalization::iSize_t e=0;e<epsilonNbr;++e,currEpsilon+=epsilonStep)
        {
            sqlRow.UpdateValue("Epsilon",currEpsilon);
            
            double currTheta = thetaMin;
        
            for(diagonalization::iSize_t t=0;t<thetaNbr;++t,currTheta+=thetaStep,++idCounter)
            {
                sqlRow.UpdateValue("Theta",currTheta);
                
                std::stringstream ss;
            
                ss.str("");

                ss<<"band_data_kx_"<<optionList["kx"].as<diagonalization::iSize_t>()<<"_ky_"<<optionList["ky"].as<diagonalization::iSize_t>()<<"_id_"<<idCounter;

                sqlRow.UpdateValue("OutFileName",_TEXT_,ss.str());
                
                data.push_back(sqlRow);

                progress.Display(idCounter-startCounter+1);
            }
        }
        
        //  Insert all data into the SQL table
        
        sql.InsertIntoTable(tableName,&data);

        sql.CopyTableToTextFile(tableName,fileName.str());
        
        utilities::cout.MainOutput()<<"\n\tDONE\n"<<std::endl;
    }
    
    return 0;
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
    const int argc, //!<	Number of characters to parse
	char *argv[])	//!<	Character array to parse
{
    namespace po = boost::program_options;

    po::variables_map vm;

    //  Construct sql options list
    po::options_description sqlOpt(diagonalization::myOptions::GetCommonSqlOptions());
    
    diagonalization::myOptions::AddGenSqlOptions(sqlOpt);

    //  Main program description
    po::options_description allOpt("\n\tThis program initializes an SQL database for the optical flux lattice model.\n\t\n\tSee e.g. PRL109,265301 (2013) for the model\n\n\tThe program input options are as follows");
        
    //	Declare option groups included
    allOpt.add(diagonalization::myOptions::GetGeneralOptions()).add(diagonalization::myOptions::GetCommonInteractingModelOptions()).add(sqlOpt);

    try
    {
        //	Map the command line arguments
        po::store(po::command_line_parser(argc,argv).options(allOpt).run(),vm);

        //	Respond to help request
        if(vm.count("help"))
        {
	        utilities::cout.MainOutput() << allOpt << "\n";
	        exit(EXIT_SUCCESS);
        }
    
        po::notify(vm);
    }
    catch(po::error& e)
    {
        utilities::cout.MainOutput()<<allOpt<<std::endl;
        
        std::cerr<<utilities::cout.HyphenLine()<<std::endl;
        
        std::cerr<<std::endl<<"\tERROR:\t"<<e.what()<<std::endl;
        
        std::cerr<<std::endl<<utilities::cout.HyphenLine()<<std::endl;
        
        exit(EXIT_SUCCESS);
    }
    
    //  Set global verbosity level

    utilities::cout.SetVerbosity(vm["verbose"].as<int>());

    utilities::cout.MainOutput()<<"\n\tRun with -h option to see program options"<<std::endl;

    return vm;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
