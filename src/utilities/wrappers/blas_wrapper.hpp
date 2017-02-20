////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 06/09/2014
//!
//!  \file
//!     A c++ wrapper around a number of selected BLAS subroutines. 
//!     All the functions are templated, with a templated-based switch 
//!     between BLAS subroutines corresponding to different variable
//!     types (principally between double and complx<double> types).
//!
//!                    Copyright (C) 2014 Simon C Davenport
//!
//!		This program is free software: you can redistribute it and/or modify
//!		it under the terms of the GNU General Public License as published by
//!		the Free Software Foundation, either version 3 of the License,
//!		or (at your option) any later version.
//!
//!		This program is distributed in the hope that it will be useful, but
//!		WITHOUT ANY WARRANTY; without even the implied warranty of 
//!		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
//!		General Public License for more details.
//!
//!		You should have received a copy of the GNU General Public License
//!		along with this program. If not, see <http://www.gnu.org/licenses/>.
//! 
////////////////////////////////////////////////////////////////////////////////

#ifndef _BLAS_WRAPPER_HPP_INCLUDED_
#define _BLAS_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../general/dcmplx_type_def.hpp"  // Define double complex type as dcmplx
#include "../general/template_tools.hpp"   // For is_same template argument checker

#if _DEBUG_
#include "../general/debug.hpp"
#endif

///////		Import selected subroutines from the BLAS library		    ////////

extern "C"
{

////////////////////////////////////////////////////////////////////////////////
//! \brief zdscal_ scales a double complex vector X with a double precision 
//! scalar ALPHA
//! 
////////////////////////////////////////////////////////////////////////////////

void zdscal_(int* N,double* ALPHA,dcmplx* X,int* INCX);

////////////////////////////////////////////////////////////////////////////////
//! \brief zdscal_ scales a real vector X with a double precision 
//! scalar ALPHA
//! 
////////////////////////////////////////////////////////////////////////////////

void dscal_(int* N,double* ALPHA,double* X,int* INCX);

////////////////////////////////////////////////////////////////////////////////
//! \brief zdotc_ is a BLAS routine for a complex double vector dot product, 
//! conjugating the first vector
//! 
////////////////////////////////////////////////////////////////////////////////

void zdotc_(dcmplx* RETURN,int* N,dcmplx* ZX,int* INCX,dcmplx* ZY,int* INCY);

////////////////////////////////////////////////////////////////////////////////
//! \brief ddot_ is a BLAS routine for a double vector dot product
//! 
////////////////////////////////////////////////////////////////////////////////

void ddot_(double* RETURN,int* N,double* DX,int* INCX,double* DY,int* INCY);

}

//  DUMMY FUNCTON DECLATIONS USED TO AVOID COMPILER ERROR WHEN USING THE
//  TEMPLATE SWITCH. THESE FUNCTIONS ARE NEVER CALLED.

void zdotc_(double* RETURN,int* N,double* ZX,int* INCX,double* ZY,int* INCY);
void ddot_(dcmplx* RETURN,int* N,dcmplx* DX,int* INCX,dcmplx* DY,int* INCY);
void ddot_(double* RETURN,int* N,dcmplx* DX,int* INCX,dcmplx* DY,int* INCY);
void zdscal_(int* N,double* ALPHA,double* X,int* INCX);
void dscal_(int* N,double* ALPHA,dcmplx* X,int* INCX);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

namespace utilities
{

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function Namespace for linear algebra routines
//!
////////////////////////////////////////////////////////////////////////////////

namespace linearAlgebra
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////////
    //!	\brief Calculate the inner product of two vectors A.B
    //!
    //! This function calls an underlying level 1 BLAS subroutine to do the 
    //! calculation for either double or complex double types.
    //!
    //! For complex variables, this function calculates the standard inner
    //! product A*.B
    //!
    //! \return the value of the inner product <A|B>
    //!
    //////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    T DotProduct(
        T* A,   //!<    Pointer to start of A vector
        T* B,   //!<    Pointer to start of B vector
        int N)  //!<    Dimension of A and B vectors
    {
        T returnValue = 0;
        int one = 1;    //  Set the vector increment value to be 1 

        if(is_same<T,dcmplx>::value)
        {
            zdotc_(&returnValue,&N,A,&one,B,&one);  
        }
        else if(is_same<T,double>::value)
        {
            ddot_(&returnValue,&N,A,&one,B,&one);
        }
        else
        {
            //  Otherwise perform a backup implementation 

            for(int i=0;i<N;i++)
            {
                returnValue += A[i]*B[i];
            }
        }
        return returnValue;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////////
    //!	\brief Scale a vector by a real scale factor
    //!
    //////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    void ScaleVector(
        T* X,           //!<    Pointer to start of vector
        double ALPHA,   //!<    Scale factor
        int N)          //!<    Vector dimension
    {
        int one = 1;    //  Set the vector increment value to be 1 

        if(is_same<T,dcmplx>::value)
        {
            zdscal_(&N,&ALPHA,X,&one);  
        }
        else if(is_same<T,double>::value)
        {
            dscal_(&N,&ALPHA,X,&one); 
        }
        else
        {
            //  Otherwise perform a backup implementation 

            for(int i=0;i<N;i++)
            {
                X[i] *= ALPHA;
            }
        }
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////////
	//!	\brief This function evaluates the mod squared overlap matrix of two lists 
	//! of vectors to give a rectangular overlap matrix of dimension dimLeft*dimRight
	//!
	//////////////////////////////////////////////////////////////////////////////////

    template<typename T>
	void OverlapMatrix(
	    double* matrix,      //!<   Output matrix (preallocated dimension dimLeft x dimRight)
	    T* leftVectors,      //!<   Array containing list of vectors back to back
	    const int dimLeft,   //!<   Number of vectors contained in the left list
	    T* rightVectors,     //!<   Array containing list of vectors back to back
	    const int dimRight,  //!<   Number of vectors contained in the right list
	    int vecDim)          //!<   Dimension of each vector
	{

	    double* p_matrix = matrix;
	    T* p_left = leftVectors;
	    
        for(int i=0;i<dimLeft;i++,p_left+=vecDim)
        {
            T* p_right  = rightVectors;
            
            for(int j=0;j<dimLeft;j++,p_right+=vecDim,p_matrix++)
            {
                //  Call level-1 BLAS to perform the dot products
            
                double temp;
                int one = 1;
                
                if(is_same<T,dcmplx>::value)
                {
                    dcmplx output;

                    zdotc_(&output,&vecDim,p_left,&one,p_right,&one);

                    temp = abs(output);
                    
                    *p_matrix = temp*temp;
                }
                else if(is_same<T,double>::value)
                {
                    double output;
                
                    ddot_(&output,&vecDim,p_left,&one,p_right,&one);
                    
                    temp = output;
                }
                
                //  Take the mod squared of each matrix element
                
                *p_matrix = temp*temp;
            }
        }
	    
	    return;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace linearAlgebra 

}   //  End namespace utilities

#endif

