////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport 
//!                                                                             
//!	 \file
//!     This file declares program options associated with observables that
//!     can be calculated/plotted for the interacting optical flux lattice
//!     model                                                     
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

#ifndef _OBSERVABLES_OPTIONS_HPP_INCLUDED_
#define _OBSERVABLES_OPTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <boost/program_options.hpp>
#include "../../utilities/general/orbital_and_state_defs.hpp"

namespace diagonalization
{
namespace myOptions
{
    namespace po = boost::program_options;
    inline po::options_description GetObservablesOptions()
    {
        po::options_description obsOpt("Observables Options");
	    obsOpt.add_options()
	    ("calculate-density-density", po::bool_switch()->default_value(0),
	     "Set to generate the density-density function and store it in a file.\n")
	    ("calculate-participation-ratio", po::bool_switch()->default_value(0),
	     "Set to generate the participation ratio (1/sum_i |psi_i|^4) and store it in a file.\n")
	    ("calculate-occupations", po::bool_switch()->default_value(0),
	     "Set to generate the probability distribution of interacting orbitals for each eigenvector.\n")
	    ("calculate-susceptibility", po::bool_switch()->default_value(0),
	     "Set to calculate the mean squared magnetization (which gives a measure of the susceptibility).\n")
	    ("most-probable-list", po::bool_switch()->default_value(false),
	     "Set to generate a list of the most probable eigenstates\n")
	    ("nbr-most-probable", po::value<iSize_t>()->default_value(1),
	     "Set the number of most probable eigenstates to store when the most-probable-list option is enabled.\n")
	    ("calculate-translational-density-density", po::bool_switch()->default_value(0),
	     "Set to calculate the translational density-density function sum_{k1,k2} <c^+_{k2-G}c_k2c^+_{k1+G}c_k1> as a function of G.\n")
        ("calculate-rotational-density-density", po::bool_switch()->default_value(0),
	     "Set to calculate the rotational density-density function sum_{k1,k2} <c^+_{~R60 k2}c_k2c^+_{R60 k1}c_k1> as a function of the 6-fold rotation.\n");
	    return obsOpt;
    };
    
    inline po::options_description GetInteractingModelPlotOptions()
    {
        po::options_description plotOpt("Plot Options");
	    plotOpt.add_options()
	    ("plot-hamiltonian", po::bool_switch()->default_value(0),
	     "Set to generate matrix plots of the interacting Hamiltonian (if it has been diagonalized).\n")
	    ("plot-occupations", po::bool_switch()->default_value(0),
	     "Set to generate a plot of the probability of orbital occupations.\n");
	    return plotOpt;
    };
}   //  End namespace myOptions
}   //  End namespace diagonalization

#endif
