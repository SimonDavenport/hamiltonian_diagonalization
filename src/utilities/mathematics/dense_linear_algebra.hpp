////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		Various dense linear algebra functions. This .h file contains a number
//!     a function templates. Several additional functions are implemented
//!     in the file linearAlgebra.cpp.
//!
//!                    Copyright (C) Simon C Davenport
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

#ifndef _DENSE_LINEAR_ALGEBRA_HPP_INCLUDED_
#define _DENSE_LINEAR_ALGEBRA_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <string.h>
#include "../general/dcmplx_type_def.hpp"
#include "../general/pi_const_def.hpp"
#include "../general/template_tools.hpp"
#if _DEBUG_
#include "../general/debug.hpp"
#endif
#if _ENABLE_HIGH_PRECISION_
#include "../wrappers/high_precision_wrapper.hpp"
#endif

namespace utilities
{
namespace linearAlgebra
{

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_HIGH_PRECISION_

	//////////////////////////////////////////////////////////////////////////////////
	//!	\brief This function performs a LU decomposition of the matrix A using the
	//!	Crout algorithm. This template uses arbitrary precision functions.
	//!
	//!	The routine takes the high precision wrapper class as an argument
	//!	P is the precision for high precision variables (r.g. mpc_t)
	//!	
	//!	Note that for simplicity of template construction, the abs value for 
	//!	complex quantities is defined as the real part of a complex number 
	//!	with 0 imaginary part
	//!		
	//!	Algorithm modified from 
	//!	http://www.mymathlib.com/matrices/linearsystems/crout.html
	//!	
	//!	Note that the input matrix A will be overwritten with the LU output.
	//!	
	//!	Template argument U corresponds to the arbitrary precision type.
	//!
	//!	Template argument P corresponds to the precision level. 
	//!	
	//!	\return The number of pivoting operators that took place
	//!	
	//////////////////////////////////////////////////////////////////////////////////

	template <typename U,int P>
	int CroutLuDecomposition(
		HpWrap<U,P>* A, 	//!<	Pointer to stored matrix (using HpWrap wrapper)
		const int dim,		//!<	Dimension of matrix	
		int* pivot)			//!<	Pivot vector: stores a list of the row	permutations
							//!		when pivoting occurs: i-th element is pivot row 
							//!		interchanged with row i	
	{

		//	Declare local variables	
		int i,j,k;							//	Generic loop variables
		HpWrap<U,P> *p_k=0,*p_row=0,*p_col=0;//	Pointers to a row k, or a given row or 
											//	column of the LU matrix
		HpWrap<U,P> myMax;					//	Max absolute value of a LU element - use
											//	to check for pivoting
		HpWrap<U,P> tmp;					//	a temporary value	
		HpWrap<U,P> myAbs,myAbs1;			//	Absolute values of elements of LU rows -
											//	use to check for pivoting
		int pivotCount=0;					//	Counts the number of pivot permutations
											//	to keep track of the sign of LU
											
		//	For each row and column, k = 0, ..., dim-1,

		for(k = 0, p_k = A; k < dim; p_k += dim, k++) {

			//	Find the pivot row

			pivot[k] = k;myMax=*(p_k + k);
			for(j = k + 1, p_row = p_k + dim; j < dim; j++, p_row += dim)
			{
				myAbs.AbsOf(myMax);
				myAbs1.AbsOf(*(p_row+k));

				if(myAbs<myAbs1)
				{
					myMax=*(p_row + k);pivot[k] = j;p_col = p_row;
				}
			}

			//	And if the pivot row differs from the current row, then interchange the two rows.

			if(pivot[k] != k){pivotCount++;
				for (j = 0; j < dim; j++) {
					myMax=*(p_k + j);
					*(p_k + j)=*(p_col + j);
					*(p_col + j)=myMax;
				}
			}

			//	Find the upper triangular matrix elements for row k.

			for(j = k+1; j < dim; j++) {*(p_k + j)/=*(p_k + k);}

			//	Update remaining matrix

			for(i = k+1, p_row = p_k + dim; i < dim; p_row += dim, i++)
			{
				for(j = k+1; j < dim; j++)
				{
					tmp=*(p_row + k);
					tmp*=*(p_k + j);
					*(p_row + j)-=tmp;;
				}
			}
		}

		return pivotCount;
	}
	
#endif
////////////////////////////////////////////////////////////////////////////////	

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	//////////////////////////////////////////////////////////////////////////////////
	//!	\brief This function performs a LU decomposition of the matrix A using the
	//!	Crout algorithm. This template uses type U variables. 
	//!
	//!	Algorithm modified from 
	//!	http://www.mymathlib.com/matrices/linearsystems/crout.html
	//!	
	//!	Note that the input matrix A will be overwritten with the LU output.
	//!	
	//!	Template argument U corresponds to the type of the variables
	//!	
	//!	\return The number of pivoting operators that took place
	//!	
	//////////////////////////////////////////////////////////////////////////////////

	template <typename U>
	int CroutLuDecomposition(
		U *A,			//!<	Type U template variable
		const int dim,	//!<	Dimension of matrix		
		int* pivot)		//!<	Pivot vector: stores a list of the row	permutations
						//!		when pivoting occurs: i-th element is pivot row 
						//!		interchanged with row i	
	{
		//	Declare local variables	
		int i,j,k;							//	Generic loop variables
		U *p_k,*p_row,*p_col;				//	Pointers to a row k, or a given row or 
											//	column of the LU matrix	
		U myMax;							//	Max absolute value of a LU element - use
											//	to check for pivoting
		double myAbs,myAbs1;				//	Absolute values of elements of LU rows -
											//	use to check for pivoting
		int pivotCount=0;					//	Counts the number of pivot permutations
											//	to keep track of the sign of LU	
										
		//	For each row and column, k = 0, ..., dim-1,

		for (k = 0, p_k = A; k < dim; p_k += dim, k++) {

			//	Find the pivot row

			pivot[k] = k;myMax=*(p_k + k);
			for (j = k + 1, p_row = p_k + dim; j < dim; j++, p_row += dim)
			{
				myAbs=abs(myMax);
				myAbs1=abs(*(p_row+k));

				if(myAbs<myAbs1)
				{
					myMax=*(p_row + k);pivot[k] = j;p_col = p_row;
				}
			}

			//	And if the pivot row differs from the current row, then interchange the two rows.

			if (pivot[k] != k){pivotCount++;
				for (j = 0; j < dim; j++) {
					myMax=*(p_k + j);
					*(p_k + j)=*(p_col + j);
					*(p_col + j)=myMax;
				}
			}

			//	Find the upper triangular matrix elements for row k.

			for (j = k+1; j < dim; j++) {*(p_k + j)/=*(p_k + k);}

			//	Update remaining matrix

			for (i = k+1, p_row = p_k + dim; i < dim; p_row += dim, i++)
			{
				for (j = k+1; j < dim; j++)
				{
					*(p_row + j)-=*(p_row + k)**(p_k + j);
				}
			}
		}

		return pivotCount;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_HIGH_PRECISION_	
	
	//////////////////////////////////////////////////////////////////////////////////
	//!	\brief Calculate the log of the determinant of a matrix stored with HpWrap class
	//!	Variables. 
	//!
	//!	This is done using LU decomposition then taking the trace.
	//!	
	//!	The routine takes the high precision wrapper class as an argument
	//!	P is the precision for high precision variables (r.g. mpc_t)
	//!		
	//!	Note that the input matrix A will be overwritten with the LU output	
	//!	
	//!	Template argument U corresponds to the arbitrary precision type.
	//!
	//!	Template argument P corresponds to the precision level. 
	//!
	//!	\return The (complex) value of the log of determinant.
	//!
	//////////////////////////////////////////////////////////////////////////////////	
		
	template <typename U,int P>
	dcmplx LogDeterminant(
		HpWrap<U,P>* A,	//!<	Memory address of the stored matrix
		const int dim)		//!<	Dimension of the matrix
	{
		//	Declare local variables

		int pivotCount; 	//  Count the number of pivot permutations of the matrix.
							//  Each permutation multiplies the determinant by -1.
		int pivot[dim];	    //	Pivoting vector
		HpWrap<U,P> det;	//	To store the value of the determinant
		dcmplx returnVal;	//	Value to return

		//	First perform an LU decomposition of the matrix A

		pivotCount=CroutLuDecomposition<U,P>(A,dim,pivot);

		//	Calculate the determinant from the trace

		det.Set(1.0);

		for(int i=0;i<dim;i++)
		{
			det*=*(A+i*dim+i);
		}

		//	Take the log

		det.LogOf(det);

		returnVal=det.Get();

		//	Include contribution from pivoting steps

		returnVal+=dcmplx(0,PI*pivotCount);

		return returnVal;
	}

#endif
////////////////////////////////////////////////////////////////////////////////
	
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	//////////////////////////////////////////////////////////////////////////////////
	//!	\brief Calculates the log of the determinant of a matrix.
	//!	
	//!	This is done using LU decomposition then taking the trace.
	//!	This is an overload of the LogDeterminant function to U types.
	//!
	//!	Note that the input matrix A will be overwritten with the LU output.
	//!
	//!	\return The (complex) value of the log of determinant
	//!
	//////////////////////////////////////////////////////////////////////////////////

	template <typename U>
	dcmplx LogDeterminant(
		U* A,			//!<	Memory address of the stored matrix
		const int dim)	//!<	Dimension of the matrix
	{
		//	Declare local variables	
		int pivotCount; 	//  Count the number of pivot permutations of the matrix.
							//  Each permutation multiplies the determinant by -1.
		int pivot[dim];	    //	Pivoting vector	
		U det;				//	To store the value of the determinant
		dcmplx returnVal;	//	Value to return 

		//	First perform an LU decomposition of the matrix A

		pivotCount=CroutLuDecomposition<U>(A,dim,pivot);

		//	Calculate the determinant from the trace

		det = 1.0;

		for(int i=0;i<dim;i++)
		{
			det*=*(A+i*dim+i);
		}

		//	Take the log

		det=log(det);

		returnVal=det;

		//	Include contribution from pivoting steps

		returnVal+=dcmplx(0,PI*pivotCount);

		return returnVal;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	//////////////////////////////////////////////////////////////////////////////////
	//!	\brief This function performs a PJP^T decomposition of the matrix A for 
	//!	calculating a Pfaffian.	
	//!
	//!	In this decomposition J is a unit antisymmetric matrix of the form
	//!	J_{i,i-(-1)^i}=-(-1)^i      for i=1,2,...,N	
	//!	
	//!	P is a lower triangular matrix	
	//!	
	//!	The algorithm is modified from ARXIV 1102.3576V2, which describes
	//!	an algorithm similar to the Crout algorithm for LU decomposition.
	//!	
	//!	Pivoting scheme is not required for this algorithm because the small
	//!	quantity gets placed in	the bottom right corner only, so you don't have 
	//!	to divide by it as in the Crout algorithm for LU decomposition. 
	//!	Also, high precision is not required.	
	//!	
	//!	Note that the input matrix A will be overwritten with the PJP output.
	//!	
	//!	Template argument U corresponds to the type.
	//!
	//////////////////////////////////////////////////////////////////////////////////

	template <typename U>
	void PjpDecomposition(
		U* A,			//!<	Memory address of the stored matrix
		const int dim)	//!<	Dimension of the matrix
	{
		//	Declare Local Variables	
		int i,j,k;				//	Generic loop variables
		U *p_k,*p_k1;			//	Pointers to the k and k+1 th row in the matrix
		U *p_row;				//	Pointer to a row in the Pfaffian matrix	
		U temp,temp2;			//	Temporary stored values in matrix (for swapping matrix rows)
		
		// For k = 1, 3,..., dim-1,
		
		for (k = 0, p_k = A; k < dim; p_k += 2*dim, k+=2) 
		{
			p_k1=p_k+dim; 
			
			// swap pairs of rows and do the division
			
			temp=-*(p_k+k+1);		
			
			for (j = k+1; j < dim; j++) 
			{
				temp2=*(p_k1+j);
				*(p_k1+j)=*(p_k+j);
				*(p_k + j)=temp2/temp;
			}
			
			// Update remaining matrix
			
			for (i = k+2, p_row = p_k + 2*dim; i < dim; p_row += dim, i++)
			{
				for (j = i+1; j < dim; j++)
				{
					*(p_row + j)+=*(p_k+ j)**(p_k1 + i)-*(p_k1 + j)**(p_k + i);
				}
			}
		}
		
		return;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	//////////////////////////////////////////////////////////////////////////////////
	//!	\brief Calculates the log of the Pfaffian of a matrix.
	//!
	//!	This is done using PJP decomposition then taking the trace.
	//!	This function applies to U type arguments.
	//!
	//!	Note that the input matrix A will be overwritten with the PJP output
	//!	
	//!	Template argument U corresponds to the type.
	//!
	//!	\return The (complex) value of the log of Pfaffian  
	//!	
	//////////////////////////////////////////////////////////////////////////////////

	template <typename U>
	dcmplx LogPfaffian(
		U* A,			//!<	Memory address of the stored matrix
		const int dim)	//!<	Dimension of the matrix
	{
		//	Declare local variables	
		int i;				//	Generic loop variables
		U *p_row;			//	Pointer to a row in the matrix	
		U pfaff;			//	To store the value of the pfaffian
		dcmplx returnVal;	//	Value to return 

		//	First Set diagonal elements to be 1,0,1,...
		//	This renders the factorization to be unique
			
		for(i=0,p_row=A;i<dim;i+=2,p_row+=2*dim)
		{
			*(p_row+i)=dcmplx(1,0);
			*(p_row+i+dim+1)=dcmplx(0,0);
		}
		
		//	Then perform an PJP decomposition of the matrix A

		PjpDecomposition<U>(A,dim);

		//	Calculate the determinant from the trace

		pfaff = 1.0;

		for(int i=0;i<dim;i++)
		{
			pfaff*=*(A+i*dim+i);
		}

		//	Take the log

		pfaff=log(pfaff);

		returnVal=pfaff;

		return returnVal;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//							
	
	////////////////////////////////////////////////////////////////////////////////
	//!	\brief This routine calculates the inverse of the lower triangular matrix L.
	//!
	//! The super-diagonal part of the matrix is not addressed.  
	//!
	//! The algorithm follows: Let M be the inverse of L, then 
	//!
	//!	L M = I, 
	//!
	//! M[i][i] = 1.0 / L[i][i] for i = 0, ..., dim-1, and
	//!
	//! M[i][j] = -[(L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j])] / L[i][i],
	//! for i = 1, ..., dim-1, j = 0, ..., i - 1.
	//!
	//!	Source http://www.mymathlib.com/matrices/linearsystems/triangular.html
	//!	
	//!	Template argument T corresponds to the type.
	//!
	//!  \return
	//!   -  0  Success 
	//!   - -1  Failure - The matrix L is singular.
	//! 
	////////////////////////////////////////////////////////////////////////////////	
		
	template <typename T>
	int LowerTriangularInverse(
		T *L,				//!<	On input, the pointer to the first element of the matrix
							//!		whose lower triangular elements form the matrix which is 
							//!		to be inverted. On output, the lower triangular part is
							//!		replaced by the inverse.  The super-diagonal elements are
							//!		not modified. 
		const int dim)		//!<	The number of rows and/or columns of the matrix L.  
	{
		int i, j, k;
		T *p_i, *p_j, *p_k;
		T sum;

		//	Invert the diagonal elements of the lower triangular matrix L.

		for (k = 0, p_k = L; k < dim; p_k += (dim + 1), k++) {
			if (*p_k == 0.0) return -1;
			else *p_k = 1.0 / *p_k;
		}

		//	Invert the remaining lower triangular matrix L row by row.

		for (i = 1, p_i = L + dim; i < dim; i++, p_i += dim) 
		{
			for (j = 0, p_j = L; j < i; p_j += dim, j++) 
			{
				sum = 0.0;
				
				for (k = j, p_k = p_j; k < i; k++, p_k += dim)
				{
					sum += *(p_i + k) * *(p_k + j);
				}
				
				*(p_i + j) = - *(p_i + i) * sum;
			}
		}
	  
		return 0;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//								
	
	////////////////////////////////////////////////////////////////////////////////
	//!	\brief This routine calculates the inverse of the upper triangular matrix U.
	//!
	//! The sub-diagonal part of the matrix is not addressed. 
	//!
	//!	The algorithm follows: Let M be the inverse of U, then 
	//!
	//!	 U M = I, 
	//!
	//! M[dim-1][dim-1] = 1.0 / U[dim-1][dim-1] and 
	//!
	//! M[i][j] = -( U[i][i+1] M[i+1][j] + ... + U[i][j] M[j][j] ) / U[i][i], 
	//! for i = dim-2, ... , 0,  j = dim-1, ..., i+1. 
	//!
	//!	Source http://www.mymathlib.com/matrices/linearsystems/triangular.html
	//!	
	//!	Template argument T corresponds to the type.
	//!
	//!  \return
	//!   -  0  Success 
	//!   - -1  Failure - The matrix U is singular.
	//!	
	////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	int UpperTriangularInverse(
		T *U,				//!<	On input, the pointer to the first element of the matrix 
							//!		whose upper triangular elements form the matrix which is 
							//!		to be inverted. On output, the upper triangular part is
							//!		replaced by the inverse. The sub-diagonal elements are
							//!		not modified. 
		const int dim)		//!<	The number of rows and/or columns of the matrix U.
	{
		
		int i, j, k;
		T *p_i, *p_j, *p_k;
		T sum;
		
		//	Invert the diagonal elements of the upper triangular matrix U.
		
		for (k = 0, p_k = U; k < dim; p_k += (dim + 1), k++) 
		{
		   if (*p_k == 0.0) return -1;
		   else *p_k = 1.0 / *p_k;
		}
		
		//	Invert the remaining upper triangular matrix U.
		
		for (i = dim - 2, p_i = U + dim * (dim - 2); i >=0; p_i -= dim, i-- ) 
		{
		  for (j = dim - 1; j > i; j--)
		  {
			 sum = 0.0;
			 
			 for (k = i + 1, p_k = p_i + dim; k <= j; p_k += dim, k++ ) 
			 {
				sum += *(p_i + k) * *(p_k + j);
			 }
			 
			 *(p_i + j) = - *(p_i + i) * sum;
		  }
		}
		
		return 0;
	}
	
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//							
	
	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Post multiply the nrows x ncols matrix A by the ncols x mcols matrix 
	//!	B to form the nrows x mcols matrix C, i.e. C = A B.
	//!
	//!	The matrix C should be declared as T C[nrows][mcols] in the
	//!	calling routine.  The memory allocated to C should not include any
	//!	memory allocated to A or B. 
	//!	
	//!	Template argument T corresponds to the type.
	//!
	//!	Source http://www.mymathlib.com/matrices/linearsystems/triangular.html
	//!	
	////////////////////////////////////////////////////////////////////////////////
		
	template <typename T>
	void MultiplyMatrices(
		T *C,				//!<	Pointer to the first element of the matrix C.
		T *A, 				//!<	Pointer to the first element of the matrix A. 
		const int nrows,	//!<	The number of rows of the matrices A and C. 
		const int ncols,	//!<	The number of columns of the matrices A and the
							//!		number of rows of the matrix B. 
		T *B,				//!<	Pointer to the first element of the matrix B.
		const int mcols)	//!<	The number of columns of the matrices B and C.
	{
		
		T *pB;
		T *p_B;
		int i,j,k;

		for (i = 0; i < nrows; A += ncols, i++)
		{
			for (p_B = B, j = 0; j < mcols; C++, p_B++, j++)
			{
				pB = p_B;
				*C = 0.0;

				for (k = 0; k < ncols; pB += mcols, k++)
				{
					*C += *(A+k) * *pB;
				}
			}
		}
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Calculates the matrix inverse of matrix M and overwrites the original
	//!	matrix M with its inverse.
	//!
	//!	To do the matrix inverse, first perform and LU decomposition of the matrix
	//!	PA=LU where P is a permutation matrix.
	//!	Then perform the matrix inverse of L and U to construct M^-1 = P^dagger U^-1 L^-1
	//!	N.B we need to undo the pivoting operations also
	//!
	//!	M must be square matrix of dimension dim.
	//!	
	//!	Template argument T corresponds to the type.
	//!
	////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	void DenseMatrixInverse(
		T *M,			//!<	Memory address of the matrix to be inverted
		const int dim)	//!<	Dimension of the matrix
	{
		int pivot[dim];	    //	Pivoting vector		
		//	The i-th element is the pivot row interchanged with row i
		
		//	First perform an LU decomposition
		
		CroutLuDecomposition<T>(M,dim,pivot);
		
		T *upperMatrix = new T[dim*dim];
		T *lowerMatrix = new T[dim*dim];
		
		int j1,j2;
		
		for(j1=0;j1<dim;j1++)
		{
			for(j2=0;j2<dim;j2++)
			{
				if(j2>j1)
				{
					upperMatrix[j1*dim+j2]=M[j1*dim+j2];
					lowerMatrix[j1*dim+j2]=0.0;
				}
				else if(j1==j2)
				{
					//upperMatrix[j1*dim+j2]=M[j1*dim+j2];
					//lowerMatrix[j1*dim+j2]=1.0;
					
					upperMatrix[j1*dim+j2]=1.0;
					lowerMatrix[j1*dim+j2]=M[j1*dim+j2];
				}
				else
				{
					lowerMatrix[j1*dim+j2]=M[j1*dim+j2];
					upperMatrix[j1*dim+j2]=0.0;
				}
			}
		}
		
		UpperTriangularInverse(upperMatrix,dim);
		LowerTriangularInverse(lowerMatrix,dim);

		MultiplyMatrices(M,upperMatrix,dim,dim,lowerMatrix,dim);
		
		//	Now we have L^-1 U^-1 but the final step is to undo the pivoting operation
		
		dcmplx exchange,*p_col;
		
		for (j1 = dim-1; j1>=0; j1--) 
		{
			if (pivot[j1] != j1)
			{
				for (j2 = 0,p_col=M; j2 < dim; j2++,p_col+=dim) 
				{
					exchange=*(p_col+j1);
					*(p_col+j1)=*(p_col+pivot[j1]);
					*(p_col+pivot[j1])=exchange;
				}
			}
		}
		
		delete[] upperMatrix;
		delete[] lowerMatrix;
		
		return;
	}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//	

    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief Replaces a rectangular complex matrix with its complex conjugate.
    //!
    //! This function only performs an operation when the tempalte parameter
    //! is a complex type.
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    void MatrixConjugate(
        dcmplx* M,		    //!<	Memory address of the matrix
        const int dim1,		//!<	Row dimension of the matrix
        const int dim2)		//!<	Column dimension of the matrix
    {
        if(is_same<T,dcmplx>::value)
        {
            dcmplx *p_M;	//	Pointer to a matrix element

            p_M=M;

            for(int i=0;i<dim1;i++)
            {
                for(int j=0;j<dim2;j++,p_M++)
                {
                    *p_M=std::conj(*p_M);
                }
            }
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//	
    
    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief Normalises a complex vector such that V^+ V = 1 or a real vector
    //! such that V^T V = 1
    //!
    //! Overwrites the original vector with the normalised version
    //! High precision variables can optionally be used to prevent double variable 
    //! overflow. High precision is also automatically used in the case of any nan
    //! or inf results. 
    //!
    //! \return The norm of the vector 
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    double NormalizeVector(
        T* V,				    //!<	Memory address of an array storing the 
                                //!		vector
        const int dim,			//!<	Dimension of the vector
        bool useHighPrec)		//!<	Option to accumulate total value with 
                                //!		high precision variables
    {
        //	Declare local variables
        T *p_V = V;
        double total=0.0;
        double norm;

        //	First try with standard precision functions

        if(is_same<T,dcmplx>::value)
        {
            for(int i=0;i<dim;i++,p_V++)
            {
                total += real( (*p_V) * std::conj(*p_V));
            }
        }
        else if(is_same<T,double>::value)
        {
            for(int i=0;i<dim;i++,p_V++)
            {
                total += *p_V * *p_V;
            }
        }
        
        norm=sqrt(total);

        ////////////////////////////////////////////////////////////////////////////////
        #if 0
        
        //	If there was an issue, use a high precision implementation

        if(std::isnan(norm) || std::isinf(norm))
        {
            useHighPrec=true;
        }

        if(useHighPrec)
        {
            int prec=128;
        
            if(is_same<T,dcmplx>::value)
            {
                mpfr_t convert;
                mpfr_t mpfrTot;

                mpfr_init2(convert,prec);
                mpfr_init2(mpfrTot,prec);

                mpfr_set_d(mpfrTot,0.0,MPFR_RNDN);

                p_V=V;

                for(int i=0;i<dim;i++,p_V++)
                {
                    mpfr_set_d(convert,real(*p_V*conj(*p_V)),MPFR_RNDN);

                    mpfr_add(mpfrTot,mpfrTot,convert,MPFR_RNDN);
                }

                mpfr_sqrt(mpfrTot,mpfrTot,MPFR_RNDN);

                norm=mpfr_get_d(mpfrTot,MPFR_RNDN);

                mpfr_clear(convert);
                mpfr_clear(mpfrTot);
            }
            else if(is_same<T,double>::value)
            {
                //  TODO
            }
        }

        #endif
        ////////////////////////////////////////////////////////////////////////////////
        
        //	Re-write the vector with the normalised values

        p_V=V;

        for(int i=0;i<dim;i++,p_V++)
        {
            *p_V/=norm;
        }

        return norm;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//    
    
    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief This function calculate the overlap between two complex vectors
    //!	leftVec and rightVec of length dim as <left|right>.	
    //!
    //!	High precision variables are optionally used and automatically used in case
    //!	of nan or inf results.
    //!
    //!	\return The value of the overlap.
    //! 
    ////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    T VectorOverlap(
        T* leftVec,		        //!<	Memory address of the left vector
        T* rightVec,		    //!<	Memory address of the right vector
        const int dim,			//!<	Dimension of the vectors
        bool useHighPrec)		//!<	Option to accumulate total value with 
                                //!		high precision variables
    {
        //	Declare local variables
        T *p_left,*p_right;
        T overlap=0.0;
        T tmp;

        //	First try with standard precision functions

        if(is_same<T,dcmplx>::value)
        {
            p_left=leftVec;
            p_right=rightVec;
            for(int i=0;i<dim;i++,p_left++,p_right++)
            {
                overlap += conj(*p_left)**p_right;
            }
        }
        else if(is_same<T,double>::value)
        {
            p_left=leftVec;
            p_right=rightVec;
            for(int i=0;i<dim;i++,p_left++,p_right++)
            {
                overlap += *p_left* *p_right;
            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        #if 0
        
        //	If there was an issue, use a high precision implementation

        if(std::isnan(abs(overlap)) || std::isinf(abs(overlap)))
        {
            useHighPrec=true;
        }

        if(useHighPrec)
        {
            int prec=128;

            if(is_same<T,dcmplx>::value)
            {
                mpfr_t convert_real;
                mpfr_t convert_imag;
                mpfr_t overlap_real;
                mpfr_t overlap_imag;

                mpfr_init2(convert_real,prec);
                mpfr_init2(convert_imag,prec);
                mpfr_init2(overlap_real,prec);
                mpfr_init2(overlap_imag,prec);

                mpfr_set_d(overlap_real,0.0,MPFR_RNDN);
                mpfr_set_d(overlap_imag,0.0,MPFR_RNDN);

                p_left=leftVec;
                p_right=rightVec;
                for(int i=0;i<dim;i++,p_left++,p_right++)
                {
                    tmp=conj(*p_left)**p_right;

                    mpfr_set_d(convert_real,real(tmp),MPFR_RNDN);
                    mpfr_set_d(convert_imag,imag(tmp),MPFR_RNDN);

                    mpfr_add(overlap_real,overlap_real,convert_real,MPFR_RNDN);
                    mpfr_add(overlap_imag,overlap_imag,convert_imag,MPFR_RNDN);
                }

                overlap=dcmplx(mpfr_get_d(overlap_real,MPFR_RNDN),mpfr_get_d(overlap_imag,MPFR_RNDN));

                mpfr_clear(convert_real);
                mpfr_clear(convert_imag);
                mpfr_clear(overlap_real);
                mpfr_clear(overlap_imag);
            }
            else if(is_same<T,double>::value)
            {
                //  TODO
            } 
        }
        
        #endif
        ////////////////////////////////////////////////////////////////////////////////

        return overlap;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief Function to orthogonalise a set of nbr vectors Vin of dimension dim
    //!
    //!	  The vectors should be pre normalised.
    //!
    //!	  I'm using a simple Gram-Schmidt procedure for the orthogonalisation.
    //!	  This is not numerically stable! To get around that, re-orthogonalise 
    //!	  at least once. Note that the orthogonalised wave functions are
    //!	  normalised at the end of this routine. The buffer Vorth will contain
    //!	  the set of orthogonalised vectors at the end
    //!
    //!	  The first orthogonal vector is identical to the first vector. |0>'=|0>
    //!	  The second orthogonal vector is given by |1>'=|1> - <0|1>|0>'	
    //!	  The values of |1>',|2>' ... are buffered and then the |1>,|2>... 
    //!	  files are overwritten at the end of the calculation.
    //!	  In general |j>' = |j>-sum _ k=1 ^ j-1 '<k|j> |k>'	
    //!	  where both the |k>' and |k> vectors are normalised.
    //!
    //!	  \return Algorithm returns 1 if the wave functions are not orthogonal to 
    //!	  with the tolerance <j|k><tol. 0 is returned otherwise.
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    int GramSchmidtOrthogonalize(
        T* Vin,			        //!<	Input array containing a list of all vectors
        T* Vorth,			    //!<	Output array containing a list of orthogonal vectors
        const int nbr,			//!<	Number of vectors in the input/output array
        const int dim,			//!<	Length of all vectors
        const double tol,		//!<	Tolerance for successful return 
        bool useHighPrec)		//!<	Option to accumulate total value with 
                                //!		high precision variables
    {
        //	Declare local variables

        T overlap;
        int flag=0;
        double norm;
        T *p_orth,*p_in;

        //	Check vectors are normalised, and if not then do so now

        for(int j=0;j<nbr;j++)
        {
            //	Quick low precision check
            norm=NormalizeVector<T>(Vin+j*dim,dim,false);

            if(norm-1.0<tol)
            {
                norm=NormalizeVector<T>(Vin+j*dim,dim,true);
            }
        }

        //	Set |j>'=|j> and initialise the buffers

        memcpy(Vorth,Vin,dim*nbr*sizeof(T));

        //	Calculate |j>' = |j>-sum _ k=1 ^ j-1 '<k|j> |k>'

        for(int j=1;j<nbr;j++)
        {
            //std::cout<<"\tvec"<<j<<std::endl;

            for(int k=0;k<j;k++)
            {
                // calculate '<k|j>

                overlap=VectorOverlap<T>(Vorth+k*dim,Vin+j*dim,dim,useHighPrec);

                for(int m=0;m<dim;m++)
                {
                    Vorth[j*dim+m] -= overlap*Vorth[k*dim+m];
                }
            }

            //	normalize |j>'

            overlap=NormalizeVector<T>(Vorth+j*dim,dim,useHighPrec);

            overlap=NormalizeVector<T>(Vorth+j*dim,dim,useHighPrec);

            //	confirm orthogonality with all previous vectors

            for(int j=0;j<nbr;j++)
            {
                for(int k=0;k<j;k++)
                {
                    overlap=VectorOverlap<T>(Vorth+k*dim,Vorth+j*dim,dim,useHighPrec);

                    if(is_same<T,dcmplx>::value)
                    {
                        if(fabs(real(overlap))>tol)	flag=1;
                    }
                    else if(is_same<T,double>::value)
                    {
                        if(fabs(overlap)>tol)	flag=1;
                    }
                }
            }
        }

        return flag;
    }

}   //  End namespace linearAlgebra

}   //  End namespace utilities

#endif

