////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file declares program options associated with the non-interacting
//!     optical flux lattice model
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

#ifndef _NONINTERACTING_OFL_MODEL_OPTIONS_HPP_INCLUDED_
#define _NONINTERACTING_OFL_MODEL_OPTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <boost/program_options.hpp>
#include "../../utilities/general/orbital_and_state_defs.hpp"

namespace diagonalization
{

namespace myOptions
{
    namespace po = boost::program_options;
    inline po::options_description GetCommonNonInteractingOflModelOptions()
    {
        po::options_description singleParticleOpt("Noninteracting optical flux lattice model options");
        singleParticleOpt.add_options()
        ("params-file",po::value<std::string>()->default_value("parameters.dat"),
         "Specify the name of the file where the single particle model data are stored\nDefault behaviour is to look for this file.\n");
        return singleParticleOpt;
    };
    
    inline void AddNonInteractingOflModelGridOptions(po::options_description& singleParticleOpt)
    {
        singleParticleOpt.add_options()
        ("calculate-spatial-wavefunctions", po::value<iSize_t>()->default_value(0),
	     "Set to x>0 to calculate the wave function psi^sigma_{kx,ky}(r) for a real space grid r with x by x points. Set x=0 to turn off this option.\n")
	    ("grid-spacing", po::value<double>()->default_value(0.25),
	     "Set the value of the real-space grid spacing i.e. the distance represented by each 1x1 square on teh grid.\n")
	    ("sptial-basis", po::value<int>()->default_value(0),
	     "Specify the state basis:\n\t 0: Finite 2D k-space grid (kx,ky) \n\t 1: Hybrid localized Wannier basis (x,ky).\n")
	    ("tol", po::value<double>()->default_value(0.000000000001),
	     "Tolerence for change of basis calcualtion.\n")
        ("calculate-bandstructure", po::bool_switch()->default_value(false),
         "Set to calculate single particle band structure (automatically set if the plot-bandstructure option is set).\n")
        ("calculate-band-width", po::bool_switch()->default_value(false),
         "Set to calculate single particle band width as a function of lattice depth from 0.1 to 5.0 (automatically set if plo-band-width option is set).\n")
        ("calculate-magnetization", po::bool_switch()->default_value(false),
         "Set to calculate single particle magnetization (automatically set if plot-magnetization option is set).\n");
    }
    
    inline po::options_description GetNonInteractingOflModelPlotOptions()
    {
        po::options_description plotOpt("Plot Options");
	    plotOpt.add_options()
	    ("plot-bandstructure", po::bool_switch()->default_value(false),
	    "Set to make a pdf plot of the single particle band structure.\n")
	    ("nbr-bands", po::value<iSize_t>()->default_value(2),
         "Specify the number of single particle bands calculated (for making band structure plots only)\n")
	    ("plot-band-width", po::bool_switch()->default_value(false),
	    "Set to make a pdf plot of single particle band width and band gap vs lattice depth (V0) from 0.1 tp 5.0.\n")
	    ("plot-magnetization", po::bool_switch()->default_value(0),
	     "Set to generate a plot of the magnetization map\n")
	    ("x-grid", po::value<iSize_t>()->default_value(4),
	    "Set the kx discretisation used for the single particle model.\n")
	    ("y-grid", po::value<iSize_t>()->default_value(4),
	    "Set the kx discretisation used for the single particle model.\n")
	    ("x-cut", po::value<iSize_t>()->default_value(6),
         "Specify the x cut-off in the Bloch coefficient labels (single particle model plots only)\n")
        ("y-cut", po::value<iSize_t>()->default_value(6),
         "Specify the y cut-off in the Bloch coefficient labels (single particle model plots only)\n")
        ("x-grid-shift", po::value<double>()->default_value(0.0),
	      "Specify the shift in the kx value when constructing the k-space grid for the plot")
	    ("y-grid-shift", po::value<double>()->default_value(0.0),
	      "Specify the shift in the ky value when constructing the k-space grid for the plot");
	    return plotOpt;
    }
}   //  End namespace myOptions
}   //  End namespace diagonalization
#endif
