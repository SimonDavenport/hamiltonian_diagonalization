////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 06/10/2014
//!                                                                             
//!	 \file
//!     This file declares program options associated with ARPACK's internal
//!     options
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

#ifndef _ARPACK_OPTIONS_HPP_INCLUDED_
#define _ARPACK_OPTIONS_HPP_INCLUDED_

//  For command line options parsing
#include <boost/program_options.hpp>

//  For parameterSize_t definition
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
	    ("shift-invert-mode",po::value<bool>()->default_value(false),
	     "[DISABLED!] Option to solve the shift-inverted eigenvalue problem (A-SIGMA*I)^-1 x = nu x, for potentially faster convergence to lowest lying eigenvalues of Ax = lx.\n")
	    ("set-shift",po::value<double>()->default_value(0.0),
	     "Set the shift value, SIGMA, to be used in shift-invert mode. Eigenvalues closest to the shift factor will converge fastest.\n\nIMPORTANT - to obtain the lowest lying eigenvalue, set the SIGMA value very close to it, or increase the number of calculated eigenvalues so that it does not get cut off. Eigenvalues closest to the SIGMA value will be included first. \n\nNOTE - it is necessary to specify the shift value in quotes e.g. '-32.1'.\n")
        ("use-initial",po::value<bool>()->default_value(false),
          "Option to read in an initial ARPACK vector from a binary file.\n")
        ("initial-file",po::value<std::string>()->default_value("initVector"),
          "Specify the name of the initial vector file (additional information and extension .bin will be added by the program).\n")
        ("store-final",po::value<bool>()->default_value(false),
          "Option to store the final ARPACK vector in a binary file.\n")
        ("final-file",po::value<std::string>()->default_value("finalVector"),
          "Specify the name of the final vector file (additional information and extension .bin will be added by the program).\n");
	     
	    return arOpt;
    };

}    //  End namespace myOptions 

}   //  End namespace diagonalization

#endif
