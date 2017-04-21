////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file declares program options associated with ARPACK's internal
//!     options
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

#ifndef _ARPACK_OPTIONS_HPP_INCLUDED_
#define _ARPACK_OPTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <boost/program_options.hpp>
#include "../utilities/general/orbital_and_state_defs.hpp"

namespace diagonalization
{
namespace myOptions
{
    namespace po = boost::program_options;
    inline boost::program_options::options_description GetArpackOptions()
    {
        po::options_description arOpt("ARPACK options");
	    arOpt.add_options()
        ("use-initial", po::bool_switch()->default_value(false),
          "Set to read in an initial ARPACK vector from a binary file.\n")
        ("initial-file", po::value<std::string>()->default_value("initVector"),
          "Specify the name of the initial vector file (additional information and extension .bin will be added by the program).\n")
        ("store-final", po::bool_switch()->default_value(false),
          "Set to store the final ARPACK vector in a binary file.\n")
        ("final-file", po::value<std::string>()->default_value("finalVector"),
          "Specify the name of the final vector file (additional information and extension .bin will be added by the program).\n");
	    return arOpt;
    };
}    //  End namespace myOptions 
}   //  End namespace diagonalization
#endif
