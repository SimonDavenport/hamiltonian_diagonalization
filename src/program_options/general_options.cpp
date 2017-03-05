////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file declares the general program options to be used by any main
//!     program.
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
#include "general_options.hpp"

namespace diagonalization
{   
    namespace myOptions
    {    
        //!
        //! Convert file format integer code to a corresponding text value
        //!
        io::fileFormat_t GetFileFormat(
            const int formatCode,
            utilities::MpiWrapper& mpi)
        {
            if(0==formatCode)
            {
                return io::_TEXT_;
            }
            else if(1==formatCode)
            {
                return io::_BINARY_;
            }
            else
            {   
                std::cerr << "\n\tERROR Unknown file format code " << formatCode <<std::endl;
                mpi.m_exitFlag = true;
                return io::_TEXT_;
            }
        }
        
        //!
        //! Convert "method" option code to a corresponding value
        //!
        diagonalizationMethod_t GetDiagonalizationMethod(
            const int methodCode,
            utilities::MpiWrapper& mpi)
        {
            if(0==methodCode)	    
            {
                 return _FULL_;
            }
            else if(1==methodCode)
            {
                return _LANCZOS_;
            }
            else 
            {
                std::cerr << "\n\tERROR Unknown method code " << methodCode << std::endl;
                mpi.m_exitFlag = true;
                return _FULL_;
            }
        }
        
        //!
        //! Convert use-hash bool code to a corresponding value
        //!
        tableFormat_t GetTermStorageType(
            const bool useHash, 
            utilities::MpiWrapper& mpi)
        {
            if(useHash)
            {
                return _HASH_;
            }
            else
            {
                return _ARRAY_;
            }
        }
    }   //  End namespace myOptions
}   //  End namespace diagonalization
