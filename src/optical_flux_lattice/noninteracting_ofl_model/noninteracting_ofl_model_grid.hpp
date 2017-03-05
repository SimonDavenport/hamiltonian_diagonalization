////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store solutions to a non-interacting 
//!     optical flux lattice model at different k-space points.
//!     See e.g. PRL 109, 265301 (2012).
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

#ifndef _NONINTERACTING_OFL_MODEL_GRID_HPP_INCLUDED_
#define _NONINTERACTING_OFL_MODEL_GRID_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "noninteracting_ofl_model.hpp"
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/wrappers/fftw_wrapper.hpp"
#include "../../utilities/general/dcmplx_type_def.hpp"
#include "../../utilities/general/i_const_def.hpp"
#include "../../utilities/general/pi_const_def.hpp"
#include "../../utilities/wrappers/blas_wrapper.hpp"
#include "../../utilities/general/run_script.hpp"
#include "../../utilities/data_structures/multi_key_hash.hpp"
#include "../../utilities/wrappers/program_options_wrapper.hpp"
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif
#include <algorithm>
#include <string.h>
#include <iomanip>

namespace diagonalization
{
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief The NonInteractingOflModelGrid class contains functions that act
    //! on arrays of NoninteractingOflModel objects. No data is contained in the
    //! class itself, and default constructors/destructors are always used. 
    //! Writing the functions in this way allows us to group them together in 
    //! a friend class of the NoninteractingOflModel class.
    //!
    //! For instance, an array of single particle energy levels at different
    //! positions in k-space would represent a single-particle band structure
    //////////////////////////////////////////////////////////////////////////////////
    class NonInteractingOflModelGrid
    {
        private:
	    kState_t CombineIndexes(const kState_t x, const kState_t y, const iSize_t dim) const;
        public:
        void CalculateBandstructure(double* lowerBand, double* secondBand, double& bandWidth,
             double& bandGap, const bool getEnergyData, double* magnetization, const bool getMagnetizationData,
             const iSize_t dimX, const iSize_t dimY, const double offsetX, const double offsetY,
             NonInteractingOflModel* ham, utilities::MpiWrapper& mpi);  
        void PlotBandstructure(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void PlotBandWidth(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void PlotMagnetization(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void CalculateSpatialWaveFunctions(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void GenerateBlochWaveFunctionTable(const iSize_t cutOff, const iSize_t dimX, const iSize_t dimY, 
                                            const double offsetX, const double offsetY, NonInteractingOflModel* blochTable,
                                            const NonInteractingOflModelData* params, utilities::MpiWrapper& mpi);
        void GenerateWannierCoefficients(const iSize_t dimX, const iSize_t dimY, const double offsetX,
                                         NonInteractingOflModel* blochTable, dcmplx* M, const utilities::MpiWrapper& mpi);
        void CalcualteSigmaZMap(const iSize_t dimX, const iSize_t dimY, NonInteractingOflModel* blochTable, double* map);   
    };
}   //  End namespace diagonalization
#endif
