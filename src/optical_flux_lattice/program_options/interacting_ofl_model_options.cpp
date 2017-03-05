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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "interacting_ofl_model_options.hpp"

namespace diagonalization
{
    namespace myOptions
    {   
        //!
        //! Convert basis code to a flag to use Wannier basis or not
        //!
        bool GetBasisType(
            const int basisCode,
            utilities::MpiWrapper& mpi)
        {
            if(0==basisCode)
            {
                return false;
            }
            else if(1==basisCode)
            {
                return true;
            }
            else
            {
                std::cerr << "\n\tERROR Unknown basis code : " << basisCode << std::endl;
                mpi.m_exitFlag = true;
                return false;
            }
        }  
    }   //  End namespace myOptions
}   //  End namespace diagonalization
