////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 08/01/2015
//!                                                                             
//!	 \file
//!     This file declares program options associated with a FQHE
//!     Haldane pseudopotential Hamiltonian in the sphere geometry
//!                                                        
//!                    Copyright (C) 2015 Simon C Davenport
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

#ifndef _PSEUDOPOTENTIAL_HAMILTONIAN_OPTIONS_HPP_INCLUDED_
#define _PSEUDOPOTENTIAL_HAMILTONIAN_OPTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

//  For command line options parsing
#include <boost/program_options.hpp>

//  For iSize_t definition
#include "../../utilities/general/orbital_and_state_defs.hpp"

namespace diagonalization
{

namespace myOptions
{
    namespace po = boost::program_options;

    inline po::options_description GetPseudopotentialHamiltonianOptions()
    {
        po::options_description modelOpt("Pseudopotential Hamiltonian options");
	    modelOpt.add_options()
	    ("nbr-orbitals,o",po::value<iSize_t>()->default_value(9),
	     "Specify the number of z-component of orbital angular momentum orbitals. Equivalent to 2S+1, where 2S is the monopole strength\n")
	    ("nbr-particles,n",po::value<iSize_t>()->default_value(3),
	     "Specify number of particles in the state\n")
	    ("nbr-levels",po::value<iSize_t>()->default_value(1),
	     "Specify number of Landau levels (1 or 2 currently available).\n")
	    ("diagonalize,d",po::value<bool>()->default_value(false),
	     "Set to 1 to generate and diagonalize the model's Hamiltonian.\n")
	    ("method",po::value<iSize_t>()->default_value(0),
	     "Specify the diagonalization method:\n\t 0 Full (LAPACK) \n\t 1 Lanczos (ARPACK)\n")
	    ("two-body-pseudopotentials",po::value<std::vector<double> >()->multitoken(),
	     "A list of 2-body pseudopotential coefficients in order of lz value. e.g. 0.0 1.0 specifies V1=1 and all other pseudopotentials 0 (default).\n")
	    ("two-body-pseudopotentials-2ll",po::value<std::vector<double> >()->multitoken(),
	     "A list of 2nd Landau level 2-body pseudopotential coefficients in order of lz value. e.g. 0.0 1.0 specifies V1=1 and all other pseudopotentials 0 (default).\n")
	    ("lz-sectors",po::value<std::vector<iSize_t> >()->multitoken(),
	     "Specify (2x) the total angular momentum sector to be diagonalized. Multiple sectors are specified as a list of integers.\n")
	    ("block-diagonalize,b",po::value<bool>()->default_value(false),
	     "Set to 1 in order to enable block diagonalization of ALL available sectors (alternatively, specify the sectors required with the --lz-sectors option).\n")
	    ("nbr-eigenvalues,e",po::value<iSize_t>()->default_value(4),
	     "Specify number of lowest lying eigenvalues of the interacting Hamiltonian that will be calculated\n")
	    ("eigenvalues-file",po::value<bool>()->default_value(false),
	     "Set to 1 to store/retrieve eigenvalues in/from a file\n")
	    ("eigenvectors-file",po::value<bool>()->default_value(false),
	     "Set to 1 to additionally store/retrieve eigenvectors in/from a file (eigenvalues also stored/retrieved if this option is set)\n")
        ("hamiltonian-to-file",po::value<bool>()->default_value(false),
         "Set to 1 to store Hamiltonian matrix in a data file\n");

        return modelOpt;
    };

}   //  End namespace myOptions

}   //  End namespace diagonalization

#endif
