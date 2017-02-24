////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file contains common utility declarations used in sparse matrix
//!     classes.
//!
//!     NOTE 1: the sparse_matrix classes must use integer type indexing
//!     in order to interface with the SPARSE BLAS library. It is assumed
//!     that the largest matrix required will not exceed 2^31 x 2^31. 
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

#ifndef _SPARSE_MATRIX_HPP_INCLUDED_
#define _SPARSE_MATRIX_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////		 
#include <cstdint> 

namespace utilities
{
    //!
    //!  Class template pre-declaration - allows MappedSparseMatrix constructor
    //!  to be called from a CrsSparseMatrixType
    //!
    template <typename T>
    class CrsSparseMatrix;
    
    //!
    //!  Class template pre-declaration - allows CrsSparseMatrixType constructor
    //!  to be called from a MappedSparseMatrix
    //!
    template <typename T>
    class MappedSparseMatrix;

    //!
    //! Define a storage method enum
    //!
    enum storageMethod_t {_SPARSE_MAPPED_,_SPARSE_CRS_,_DENSE_};

}   //  End namespace utilities
#endif
