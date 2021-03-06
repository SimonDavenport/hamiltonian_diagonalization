////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file declares program options associated with the interacting 
//!     optical flux model Hamiltonian
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

#ifndef _INTERACTING_OFL_MODEL_OPTIONS_HPP_INCLUDED_
#define  _INTERACTING_OFL_MODEL_OPTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <boost/program_options.hpp>
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../program_options/general_options.hpp"
#include "../../utilities/wrappers/mpi_wrapper.hpp"

namespace diagonalization
{
    namespace myOptions
    {
        namespace po = boost::program_options;
        inline po::options_description GetCommonInteractingModelOptions()
        {
            po::options_description modelOpt("Interacting optical flux lattice model options");
	        modelOpt.add_options()
	        ("kx,x",po::value<iSize_t>()->default_value(2),
	         "Specify k-space discretization in the x direction\n")
	        ("ky,y",po::value<iSize_t>()->default_value(2),
	         "Specify k-space discretization in the y direction\n")
	        ("nbr,n",po::value<iSize_t>()->default_value(2),
	         "Specify number of particles in the state\n")
	        ("kx-shift",po::value<double>()->default_value(0.0),
	          "Specify the shift in the kx value when constructing the k-space grid")
	        ("ky-shift",po::value<double>()->default_value(0.0),
	          "Specify the shift in the ky value when constructing the k-space grid");
            return modelOpt;
        }
        
        inline void AddTermOptions(po::options_description& modelOpt)
        {
            modelOpt.add_options()
	        ("tol-vkkkk", po::value<double>()->default_value(0.000000000001),
	         "Set the tolerance for calculation of vkkkk matrix elements.\n")
	        ("min-k-cut", po::value<iSize_t>()->default_value(5),
	         "Set the starting value for the kx and ky cut-off used to calculate vkkkk.\n")
	        ("max-k-cut", po::value<iSize_t>()->default_value(18),
	         "Set the maximum allowed value for the kx cut-off used to calculate vkkkk. If the matrix element values have not converged to within the set tolerance by this point then the program will proceed with the current values of vkkkk.\n")
	        ("k-cut-step", po::value<iSize_t>()->default_value(2),
	         "Set the step increase in k-cut at which to test for convergence within the set tolerance compared with the previous k-cut value.\n");
        }
        
        inline void AddMoreInteractingModelOptions(po::options_description& modelOpt)
        {
            modelOpt.add_options()
	        ("basis", po::value<int>()->default_value(0),
	         "Specify the state basis:\n\t 0: Finite 2D k-space grid (kx,ky) \n\t 1: Hybrid localized Wannier basis (x,ky).\n")
	        ("sectors,s", po::value<std::vector<iSize_t> >()->multitoken(),
	         "Specify the total linear momentum sector as a list of integer pairs kx_tot ky_tot. e.g. -s 1 1 2 0 3 2 would diagonalize 3 sectors in succession. Note that for the Wannier basis, only ky_tot data is used but you still need to specify pairs of sectors.\n")
	        ("use-hash", po::value<bool>()->default_value(false),
	         "Set to 1 in order to use multi-key hash tables to store Hamiltonian coefficient tables Vkkkk, Ekk and momentum conserving lists. This option is automatically selected for the Wannier basis case.\n");
        }
        
        bool GetBasisType(const int basisCode, utilities::MpiWrapper& mpi);
    }   //  End namespace myOptions
}   //  End namespace diagonalization
#endif
