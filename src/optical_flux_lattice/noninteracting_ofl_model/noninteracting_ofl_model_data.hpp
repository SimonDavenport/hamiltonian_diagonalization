////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store data describing the non-interacting
//!     optical flux lattice model. See e.g. PRL 109, 265301 (2012).	 
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

#ifndef _NONINTERACTING_OFL_MODEL_DATA_HPP_INCLUDED_
#define _NONINTERACTING_OFL_MODEL_DATA_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/wrappers/sqlite_wrapper.hpp"
#include "../../utilities/wrappers/io_wrapper.hpp"
#include "../../utilities/wrappers/program_options_wrapper.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>
#if _DEBUG_
#include "../../general/debug.hpp"
#endif

namespace diagonalization
{
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A data structure to contain all of the single particle model
    //! parameters
    //! 
    //! See e.g. PRL 109, 265301 (2012) for the full definition of these parameters
    //////////////////////////////////////////////////////////////////////////////////
    struct NonInteractingOflModelData
    {
        double m_theta;       //!<    Polarization angle of incident beams
                              //!
        double m_V0;          //!<    Lattice depth parameter
                              //!
        double m_epsilon;     //!<    Ratio of Raman coupling and scalar potential
                              //!     (couples the different spin species in the model)
        double m_kappa;       //!<    Wave vector magnitude of incident beams
                              //!
        double m_mass;        //!<    mass of the cold atom
                              //!
        NonInteractingOflModelData();
        NonInteractingOflModelData(const NonInteractingOflModelData& other);
        NonInteractingOflModelData& operator=( const NonInteractingOflModelData& other);
        NonInteractingOflModelData(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        void ReadFromFile(const std::string fileName, utilities::MpiWrapper& mpi);
        void GenerateParametersFile() const;
        void ReadFromSql(const std::string tableName, const std::string fileName, const iSize_t sqlId, utilities::MpiWrapper& mpi);
        void MpiSynchronize(const iSize_t nodeId, const utilities::MpiWrapper& mpi);
    };
}   //  End namespace diagonalization
#endif
