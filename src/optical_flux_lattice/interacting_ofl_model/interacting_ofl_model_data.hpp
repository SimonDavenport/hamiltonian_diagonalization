////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store parameters defining an optical flux 
//!     lattice interacting Hamiltonian. See e.g. PRL 109, 265301 (2012)
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

#ifndef _INTERACTING_OFL_MODEL_DATA_HPP_INCLUDED_
#define _INTERACTING_OFL_MODEL_DATA_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/wrappers/io_wrapper.hpp"
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/wrappers/sqlite_wrapper.hpp"
#include "../../utilities/general/pi_const_def.hpp"
#include "../../utilities/wrappers/program_options_wrapper.hpp"
#include "../../program_options/general_options.hpp"
#include "../program_options/interacting_ofl_model_options.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    struct InteractingOflModelData
    {
        iSize_t m_dimX;                 //!<    Number of k-space samples in x direction
	    iSize_t m_dimY;                 //!<    Number of k-space samples in y direction
        iSize_t m_nbrParticles;         //!<    Number of particles in the system
	    iSize_t m_kxTot;                //!<    kx total linear momentum sector
	    iSize_t m_kyTot;                //!<    ky total linear momentum sector
	    double m_offsetX;               //!<    kx offset to k-space grid
	    double m_offsetY;               //!<    ky offset to k-space grid
        std::string m_outPath;          //!<    Path to output data directory
	    std::string m_inPath;           //!<    Path to input data directory
	    std::string m_outFileName;      //!<    Name of file where output data are stored
	    iSize_t m_nbrEigenvaluesToFind; //!<    A number to specify how many of the 
	                                    //!     lowest eigenvalues to find
	    double m_interactionStrength;   //!<    Interaction strength parameter
	    bool m_blockDiagonalize;        //!<    Set to true to diagonalize the Hamiltonian in 
                                        //!<    separate block-diagonal linear momentum sectors
        bool m_useSql;                  //!<    Flag to use parameters stored in an SQL database
        iSize_t m_sqlId;                //!<    Identifier for parameters in an sql database
        std::string m_sqlName;          //!<    Name of SQL file containing model data
        std::string m_sqlTableName;     //!<    Name of SQL table containing model data    
        std::string m_initialVectorFile;//!<    Name of ARPACK initial vector file
        std::string m_finalVectorFile;  //!<    Name of ARPACK final vector file
        tableFormat_t m_tableFormat;    //!<    A flag to keep track of the current table storage format
        tableFormat_t m_setTableFormat; //!<    Set the table format for storing Vkkkk, Ekk and
                                        //!     momentum conserving tables                   
        diagonalizationMethod_t m_method;//!<   Store the diagonalization method
	    bool m_useWannierBasis;         //!<    Option to use maximally localized Wannier basis
        bool m_termTablesBuilt;         //!<    Set to true once the term tables are built
        bool m_magnetizationCalculated; //!<    Flag to specify whether magnetisation map has been stored
        bool m_fockBasisBuilt;          //!<    Set to true once the Fock basis has been constructed
        bool m_hamiltonianBuilt;        //!<    Set to true once the Hamiltonian is built
        bool m_hamiltonianDiagonalized; //!<    Set to true once Hamiltonian is diagonalized
        InteractingOflModelData();
        ~InteractingOflModelData();
        InteractingOflModelData(const InteractingOflModelData& other);
        InteractingOflModelData(boost::program_options::variables_map* optionList, utilities::MpiWrapper& mpi);
        InteractingOflModelData& operator=(const InteractingOflModelData& other);
        void MpiSynchronize(const int nodeId, const utilities::MpiWrapper& mpi);
        void ReadFromFile(const std::string fileName, utilities::MpiWrapper& mpi);
        void ReadFromSql(const std::string tableName, const std::string fileName, const iSize_t sqlId, 
                         const bool diagonalize, utilities::MpiWrapper& mpi);
        void RescaleInteractionStrength(const double mass);
        std::string MakeBaseFileName(const io::io_t io) const;
        void UpdateSqlStatus(const utilities::MpiWrapper& mpi) const;
    };
}   //  End namespace diagonalization
#endif
