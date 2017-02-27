////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file declares program options associated with SQLite 
//!     functionality
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

#ifndef _SQL_OPTIONS_HPP_INCLUDED_
#define _SQL_OPTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <boost/program_options.hpp>
#include "../../utilities/general/orbital_and_state_defs.hpp"

namespace diagonalization
{
namespace myOptions
{
    namespace po = boost::program_options;
    inline po::options_description GetCommonSqlOptions()
    {
        po::options_description sqlOpt("SQLite options");
        sqlOpt.add_options()
        ("use-sql",po::value<bool>()->default_value(false),
         "Option to use model parameters stored in an sqllite file\n")
        ("sql-name",po::value<std::string>()->default_value("resultFileKey"),
         "Specify the name of an sqllite file name where the model data and output file names are stored\n")
        ("sql-table-name",po::value<std::string>()->default_value("PARAMETERS"),
         "Specify the name of an sqllite table where data are stored\n")
        ("sql-id",po::value<iSize_t>()->default_value(1),
         "Specify the id of the parameter set extracted form the sqllite file.\n");
         return sqlOpt;
    };
    
    inline void AddGenSqlOptions(po::options_description& sqlOpt)
    {
        sqlOpt.add_options()
        ("build-sql-table",po::value<bool>()->default_value(false),
         "An option to construct a SQL table containing a cut through parameter space\n")
        ("build-sql-table-offset",po::value<bool>()->default_value(false),
         "An option to construct a SQL table containing a list of different offsets for a fixed set of system parameters\n")
        ("build-sql-table-single-particle",po::value<bool>()->default_value(false),
         "An option to construct a SQL table containing a list of different single-particle model parameters, to produce results related to the single-particle model only.\n")
        ("sql-v0-min",po::value<double>()->default_value(0.0),
         "Specify the minimum V0 parameter value to put into the SQL table.\n")
        ("sql-v0-step",po::value<double>()->default_value(0.2),
         "Specify the V0 parameter step size to put into the SQL table.\n")
        ("sql-v0-nbr",po::value<iSize_t>()->default_value(6),
         "Specify the number of V0 parameter values to put in the SQL table (starting from sql-v0-min and increasing in steps of sql-v0-step.\n")
        ("sql-g-min",po::value<double>()->default_value(0.0),
         "Specify the minimum g parameter value to put into the SQL table.\n")
        ("sql-g-step",po::value<double>()->default_value(0.2),
         "Specify the g parameter step size to put into the SQL table.\n")
        ("sql-g-nbr",po::value<iSize_t>()->default_value(6),
         "Specify the number of g parameter values to put in the SQL table (starting from sql-g-min and increasing in steps of sql-g-step.\n")
        ("sql-kx-shift-min",po::value<double>()->default_value(0.0),
         "Specify the minimum value of the kx offset to put added to the SQL table.\n")
        ("sql-kx-shift-step",po::value<double>()->default_value(0.2),
         "Specify the kx offset step size to put added to the SQL table.\n")
        ("sql-kx-shift-nbr",po::value<iSize_t>()->default_value(6),
         "Specify the number of kx offset values to put in the SQL table (starting from sql-kx-shift-min and increasing in steps of sql-kx-shift-step.\n")
        ("sql-ky-shift-min",po::value<double>()->default_value(0.0),
         "Specify the minimum value of the ky offset to put added to the SQL table.\n")
        ("sql-ky-shift-step",po::value<double>()->default_value(0.2),
         "Specify the ky offset step size to put added to the SQL table.\n")
        ("sql-ky-shift-nbr",po::value<iSize_t>()->default_value(6),
         "Specify the number of ky offset values to put in the SQL table (starting from sql-ky-shift-min and increasing in steps of sql-ky-shift-step.\n")
        ("sql-completed",po::value<bool>()->default_value(false),
         "Generate a SQLite table with all flags set to a completed calculation. Used for rebuilding corrupted tables. \n");
    };
    
    inline void AddRunSqlOptions(po::options_description& sqlOpt)
    {
	    sqlOpt.add_options()
        ("sort-vectors",po::value<iSize_t>()->default_value(0),
         "Option to keep a running check for energy level cross overs as a function of model parameters, and to relabel the eigenvalues to keep the labelling consistent. The option value sets the number of eigenvalues to keep for the running check.\n");
    };
}   //  End namespace myOptions
}   //  End namespace diagonalization 

#endif
