////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 06/11/2014
//!                                                                             
//!	 \file
//!     This file defines a class to store an array of optical flux lattice 
//!     single particle Hamiltonians e.g. at different k-space points.
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

#ifndef _SINGLE_PARTILCE_HAMILTONIAN_ARRAY_HPP_INCLUDED_
#define _SINGLE_PARTILCE_HAMILTONIAN_ARRAY_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../../utilities/general/orbital_and_state_defs.hpp"
//  For iSize_t definition

//	Class defining a single particle Hamiltonian at a particular k-point
#include "single_particle_hamiltonian.hpp"

//  Include various utility functions
#include "../../utilities/wrappers/mpi_wrapper.hpp"
//	MPI functions and variables for parallelization

#include "../../utilities/wrappers/fftw_wrapper.hpp"
//  For 1D Fourier transform functions used in the change of basis 

#include "../../utilities/general/dcmplx_type_def.hpp"
//  Definition of a complex double type

#include "../../utilities/general/i_const_def.hpp"
//  Definition of a sqrt(-1) place holder

#include "../../utilities/general/pi_const_def.hpp"
//  Definition of a PI place holder

#include "../../utilities/wrappers/blas_wrapper.hpp"
//  A wrapper for BLAS "DotProduct" function

#include "../../utilities/general/run_script.hpp"
//  For executing command line scripts

#include "../../utilities/data_structures/multi_key_hash.hpp"
//  For use in storing the change of basis

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

#include <algorithm>    //  for max_element and min_element functions
#include <string.h>     //  for memcpy function
#include <iomanip>      //  for setw function

namespace diagonalization
{

//////////////////////////////////////////////////////////////////////////////////
//! \brief The SingleParticleHamiltonianArray class contains functions that act
//! on arrays of SingleParticleHamiltonian objects. No data is contained in the
//! class itself, and default constructors/destructors are always used. 
//! Writing the functions in this way allows us to group them together in 
//! a friend class of the SingleParticleHamiltonian class.
//!
//! For instance, an array of single particle energy levels at different
//! positions in k-space would represent a single-particle band structure
//!
//////////////////////////////////////////////////////////////////////////////////

class SingleParticleHamiltonianArray
{
    private:

	kState_t CombineIndexes(const kState_t x,const kState_t y,const iSize_t dim) const;

    public:

    //  A function to map out the band structure of the lowest two bands for a grid
    //  covering the Brillouin zone 
    void CalculateBandstructure(double* lowerBand,double* secondBand,double& bandWidth,
         double& bandGap,const bool getEnergyData,double* magnetization,const bool getMagnetizationData,
         const iSize_t dimX,const iSize_t dimY,const double offsetX,const double offsetY,
         SingleParticleHamiltonian* ham,utilities::MpiWrapper& mpi);

    //  A function to plot the band structure on a grid covering the first Brillouin zone     
    void PlotBandstructure(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);

    //  A function to plot the band width/band gap as a function of the model parameters
    void PlotBandWidth(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);

    //  A function to plot the magnetization pattern
    void PlotMagnetization(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);

    //  A function to calculate the real-space orbitals
    void CalculateSpatialWaveFunctions(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);

    //  A function to generate a table of Bloch wave function values for a grid
    //  of k-space points
    void GenerateBlochWaveFunctionTable(const iSize_t cutOff,const iSize_t dimX,const iSize_t dimY,const double offsetX,const double offsetY,
         SingleParticleHamiltonian* blochTable,const SingleParticleParameters* params,
         utilities::MpiWrapper& mpi);

    //  A function to generate the coefficients and phase factors required
    //  for the change the the hybrid localized Wannier basis
    void GenerateWannierCoefficients(const iSize_t dimX,const iSize_t dimY,const double offsetX,
         SingleParticleHamiltonian* blochTable,dcmplx* M,const utilities::MpiWrapper& mpi);

    //  A function to obtain a map of the sigma^z values associated with each
    //  k-space point
    void CalcualteSigmaZMap(const iSize_t dimX,const iSize_t dimY,
         SingleParticleHamiltonian* blochTable,double* map);   
};

}   //  End namespace diagonalization

#endif

