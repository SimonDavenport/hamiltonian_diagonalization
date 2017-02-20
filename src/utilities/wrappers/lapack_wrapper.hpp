////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 06/09/2014
//!
//!  \file
//!     A c++ wrapper around a number of selected LAPACK subroutines. 
//!     All the functions are templated, with a templated-based switch 
//!     between LAPACK subroutines corresponding to different variables 
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

#ifndef _LAPACK_WRAPPER_HPP_INCLUDED_
#define _LAPACK_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <iostream>                      // For std::cerr
#include "../general/dcmplx_type_def.hpp"// Define double complex type as dcmplx
#include "../general/template_tools.hpp" // For utilities::is_same template argument checker

#if _DEBUG_
#include "../general/debug.hpp"
#endif

///////		Import selected LAPACK subroutines      ////////////////////////////

extern "C"
{

////////////////////////////////////////////////////////////////////////////////
//! \brief zgeqrf_ computes a QR factorization of a complex M-by-N matrix A
//! 
//! See here for more details:
//!
//! http://www.netlib.org/lapack/explore-html/d8/d32/zgeqrf_8f.html
//! 
////////////////////////////////////////////////////////////////////////////////

void zgeqrf_(int* M,int* N,dcmplx* A,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,
             int* INFO);
             
////////////////////////////////////////////////////////////////////////////////
//! \brief zungqr_ converts the output of zgeqrf_ into the full Q matrix.          
//!
//! See here for more details:
//!
//! http://www.netlib.org/lapack/explore-html/df/d10/zungqr_8f.html
//!
////////////////////////////////////////////////////////////////////////////////
             
void zungqr_(int* M,int* N,int* K,dcmplx* A,int* LDA,dcmplx* TAU,dcmplx* WORK,
             int* LWORK,int* INFO);

////////////////////////////////////////////////////////////////////////////////
//! \brief dgeqrf_ computes a QR factorization of a real M-by-N matrix A
//! and 
//!
//! See here for more details:
//!
//! http://www.netlib.org/lapack/explore-html/d3/d69/dgeqrf_8f.html
//! 
////////////////////////////////////////////////////////////////////////////////

void dgeqrf_(int* M,int* N,double* A,int* LDA,double* TAU,double* WORK,int* LWORK,
             int* INFO);
             
////////////////////////////////////////////////////////////////////////////////
//! \brief dorgqr_ converts the output of dgeqrf_ into the full Q matrix. 
//!
//! See here for more details:
//!
//! http://www.netlib.org/lapack/explore-html/d9/d1d/dorgqr_8f.html    
//!
////////////////////////////////////////////////////////////////////////////////       
             
void dorgqr_(int* M,int* N,int* K,double* A,int* LDA,double* TAU,double* WORK,
             int* LWORK,int* INFO);
             
////////////////////////////////////////////////////////////////////////////////
//! \brief zgesvd_ is a LAPACK routine to compute the singular value 
//! decomposition (SVD) of a complex rectangular matrix. 
//! 
//! See here for more details:
//!
//! http://www.netlib.org/lapack/explore-html/d6/d42/zgesvd_8f.html
//! 
////////////////////////////////////////////////////////////////////////////////

void zgesvd_(char* JOBU, char* JOBVT, int* M, int* N, dcmplx* A,int* LDA, double* S, 
	        dcmplx* U, int* LDU, dcmplx* VT, int* LDVT,dcmplx* WORK, int* LWORK, 
			double* RWORK, int* INFO);
			
////////////////////////////////////////////////////////////////////////////////
//! \brief dgesvd_ is a LAPACK routine to compute the singular value 
//! decomposition (SVD) of a real rectangular matrix. 
//! 
//! See here for more details:
//!
//! http://www.netlib.org/lapack/explore-html/d8/d2d/dgesvd_8f.html
//! 
////////////////////////////////////////////////////////////////////////////////

void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A,int* LDA, double* S, 
	        double* U, int* LDU, double* VT, int* LDVT,double* WORK, int* LWORK, int* INFO);

//////////////////////////////////////////////////////////////////////////////////
//!  \brief dsyevr_ computes selected eigenvalues and, optionally, eigenvectors
//!  of a real symmetric matrix A.
//!	
//!  See here for more details:
//!	
//!	 http://www.netlib.org/lapack/explore-html/df/d3b/dsyevr_8f.html
//!
//////////////////////////////////////////////////////////////////////////////////

void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA,double* VL,
		double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W,double* Z,int* LDZ,
		int* ISUPPZ,double* WORK, int* LWORK,int* IWORK,int* LIWORK,int* INFO);

////////////////////////////////////////////////////////////////////////////////
//! \brief zheevr_ computes selected eigenvalues and, optionally, eigenvectors 
//!  of a complex Hermitian matrix A. 
//!
//! See here for more detials:
//!
//!	http://www.netlib.org/lapack/explore-html/d6/dee/zheevr_8f.html
//!
////////////////////////////////////////////////////////////////////////////////

void zheevr_(char* JOBZ,char* RANGE,char* UPLO,int* N,dcmplx* A,int* LDA,double* VL,
    double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W,dcmplx* Z, int* LDZ,
	int* ISUPPZ,dcmplx* WORK, int* LWORK,double* RWORK, int* LRWORK,int* IWORK,int* LIWORK,
	int* INFO);

////////////////////////////////////////////////////////////////////////////////
//! \brief dgeev_ computes for an N-by-N real nonsymmetric matrix A, the
//! eigenvalues and, optionally, the left and/or right eigenvectors.
//! 
//!  The right eigenvector v(j) of A satisfies
//!                   A * v(j) = lambda(j) * v(j)
//!  where lambda(j) is its eigenvalue.
//!  The left eigenvector u(j) of A satisfies
//!                u(j)**H * A = lambda(j) * u(j)**H
//!  where u(j)**H denotes the conjugate-transpose of u(j).
//! 
//!  The computed eigenvectors are normalized to have Euclidean norm
//!  equal to 1 and largest component real.
//!
//!	http://www.netlib.org/lapack/explore-html/d9/d28/dgeev_8f.html
//!
////////////////////////////////////////////////////////////////////////////////

void dgeev_(char* JOBVL,char* JOBVR,int* N,double* A,int* LDA,double* WR,double* WL,
    double* VL,int* LDVL,double* VR,int* LDVR,double* WORK,int* LWORK,int* INFO);
    
////////////////////////////////////////////////////////////////////////////////
//! \brief zgeev_ computes for an N-by-N complex nonsymmetric matrix A, the
//! eigenvalues and, optionally, the left and/or right eigenvectors.
//! 
//!  The right eigenvector v(j) of A satisfies
//!                   A * v(j) = lambda(j) * v(j)
//!  where lambda(j) is its eigenvalue.
//!  The left eigenvector u(j) of A satisfies
//!                u(j)**H * A = lambda(j) * u(j)**H
//!  where u(j)**H denotes the conjugate-transpose of u(j).
//! 
//!  The computed eigenvectors are normalized to have Euclidean norm
//!  equal to 1 and largest component real.
//!
//!	http://www.netlib.org/lapack/explore-html/d9/d28/dgeev_8f.html
//!
////////////////////////////////////////////////////////////////////////////////

void zgeev_(char* JOBVL,char* JOBVR,int* N,dcmplx* A,int* LDA,dcmplx* W,
    dcmplx* VL,int* LDVL,dcmplx* VR,int* LDVR,dcmplx* WORK,int* LWORK,int* INFO);


////////////////////////////////////////////////////////////////////////////////
//! \brief dgels_ - solve overdetermined or underdetermined real linear
//! systems involving an M-by-N matrix A, or its transpose,
//! using a QR or LQ factorization of A.
//!
//! http://www.netlib.org/lapack/explore-html/d8/dde/dgels_8f.html
//!
////////////////////////////////////////////////////////////////////////////////

void dgels_(char* TRANS,int* M,int* N,int* NRHS,double* A,int* LDA,
    double* B,int* LDB,double* WORK,int* LWORK,int *INFO);

}   //  End extern "C"

//  DUMMY FUNCTON DECLATIONS USED TO AVOID COMPILER ERROR WHEN USING THE
//  TEMPLATE SWITCH. THESE FUNCTIONS ARE NEVER CALLED.

void zgeqrf_(int* M,int* N,double* A,int* LDA,double* TAU,double* WORK,int* LWORK,
             int* INFO);
void dgeqrf_(int* M,int* N,dcmplx* A,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,
             int* INFO);             
void dorgqr_(int* M,int* N,int* K,dcmplx* A,int* LDA,dcmplx* TAU,dcmplx* WORK,
             int* LWORK,int* INFO);
void zungqr_(int* M,int* N,int* K,double* A,int* LDA,double* TAU,double* WORK,
             int* LWORK,int* INFO);
void zgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A,int* LDA, double* S, 
	        double* U, int* LDU, double* VT, int* LDVT,double* WORK, int* LWORK, 
			double* RWORK, int* INFO);
void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, dcmplx* A,int* LDA, double* S, 
	        dcmplx* U, int* LDU, dcmplx* VT, int* LDVT,dcmplx* WORK, int* LWORK,
	        int* INFO);
void zheevr_(char* JOBZ,char* RANGE,char* UPLO,int* N,double* A,int* LDA,double* VL,
    double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W,double* Z, int* LDZ,
	int* ISUPPZ,double* WORK, int* LWORK,double* RWORK, int* LRWORK,int* IWORK,int* LIWORK,
	int* INFO);	              
void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, int* N, dcmplx* A, int* LDA,double* VL,
		double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W,dcmplx* Z,int* LDZ,
		int* ISUPPZ,dcmplx* WORK, int* LWORK,int* IWORK,int* LIWORK,int* INFO);
void zgeev_(char* JOBVL,char* JOBVR,int* N,double* A,int* LDA,dcmplx* W,
    double* VL,int* LDVL,double* VR,int* LDVR,double* WORK,int* LWORK,int* INFO);
void dgeev_(char* JOBVL,char* JOBVR,int* N,dcmplx* A,int* LDA,double* WR,double* WL,
    dcmplx* VL,int* LDVL,dcmplx* VR,int* LDVR,dcmplx* WORK,int* LWORK,int* INFO);

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

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function calls a LAPACK routine to calculate the QR decomposition
    //! of a rectangular matrix
    //!
    //! NOTE: the matrix A is overwritten with the R matrix
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool QrDecomposition(
        T* A,                   //!<    Memory address of matrix for which the QR decomposition is required
        T* Q,                   //!<    Memory address of matrix Q to be returned
        const int nbrRows,      //!<    Number of rows in the matrix
        const int nbrCols)      //!<    Number of columns in the matrix
    {
        //////      SET PARAMETER VALUES FOR zgeqrf_/zungqr_ LAPACK ROUTINE     ////////
    
        int M      = nbrRows;   	//  Number of rows in the input matrix
        int N      = nbrCols;   	//  Number of columns in the input matrix
        int K      = std::min(M,N); //  Number of elementary reflectors used
                                    //  in the intermediate solution
        int LDA    = std::max(1,M); //  The leading dimension of the input matrix
        T* TAU = new T[std::min(M,N)];  
                                    //  Array to store elementary reflectors
        int LWORK  = std::max(1,2*N);
                                    //  Working memory space size
        T* WORK = new T[LWORK];     //  Working memory space   
        int INFO;               	//  Diagnostic return value

        if(utilities::is_same<T,dcmplx>::value)
        {
            zgeqrf_(&M,&N,A,&LDA,TAU,WORK,&LWORK,&INFO);
        }
        else if(utilities::is_same<T,double>::value)
        {
            dgeqrf_(&M,&N,A,&LDA,TAU,WORK,&LWORK,&INFO);
        }
        
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zgeqrf_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dgeqrf_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }
        
        //  Make a copy of the returned A matrix to be used to generate
        //  the Q matrix
        
        memcpy(Q,A,N*M*sizeof(T));
        
        //  Set the lower triangular entries of A to be zero so that
        //  we preserve only the R matrix
        
        for(int i=0;i<M;i++)
        {
            for(int j=0;j<i;j++)
            {
                A[j*M+i] = 0.0;
            }
        }
        
        //  Extract the final Q matrix from the currently stored Q
        
        if(utilities::is_same<T,dcmplx>::value)
        {
            zungqr_(&M,&N,&K,Q,&LDA,TAU,WORK,&LWORK,&INFO);
        }
        else if(utilities::is_same<T,double>::value)
        {            
            dorgqr_(&M,&N,&K,Q,&LDA,TAU,WORK,&LWORK,&INFO);
        }
        
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zungqr_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dorgqr_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }
        
        return false;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function calls a LAPACK routine to evaluate the singular
    //! values and singular vectors of a rectangular matrix.
    //!
    //! NOTE: the matrix A is overwritten with the right singular vectors
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool LeftSingularValueDecomposition(
        T* A,                   //!<    Memory address of matrix for which the SVD is required
        T* U,                   //!<    Memory address to store the left singular vectors 
        double* S,              //!<    Memory address of singular values
        const int nbrRows,      //!<    Number of rows in the matrix
        const int nbrCols)      //!<    Number of columns in the matrix
    {
        //////      SET PARAMETER VALUES FOR zgesvd_/dgesvd_ LAPACK ROUTINE     ////////
        
        char JOBU  = 'S';       	//  Setting to store the left singular vectors
                                    //  in the U array
        char JOBVT = 'O';       	//  Setting to store the right singular vectors
                                    //  OVERWRITE the input matrix A
        int M      = nbrRows;   	//  Number of rows in the input matrix
        int N      = nbrCols;   	//  Number of columns in the input matrix
        int LDA    = std::max(1,M); //  The leading dimension of the input matrix
        int LDU    = M;             //  Leading dimension of the array U
        T* VT = 0;                  //  Not addressed with current JOBVT setting
        int LDVT   = N;             //  Leading dimension of the array VT 
        int LWORK  = std::max(1,2*N+2*M);
                                    //  Working memory space size
        T *WORK = new T[LWORK];     //  Working memory space
        int INFO;               	//  Diagnostic return value

        //  Change from column major to row major indexing
        //  in order to correctly pass to fortran routine	
        
        //dcmplx* B  = new dcmplx[nbrRows*nbrCols];//   Matrix to calculate the SVD for
        
        //for(int i=0;i<nbrRows;i++)
        //{
        //    for(int j=0;j<nbrCols;j++)
        //    {   
        //        B[i*nbrCols+j] = A[j*nbrRows+i];
        //    }
        //}

        if(utilities::is_same<T,dcmplx>::value)
        {
            double *RWORK = new double[LWORK];
                                    //  Working memory space
        
            //////      CALL LAPACK zgesvd_ ROUTINE     ////////////////////////////
        	
            zgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,
                &LDVT,WORK,&LWORK,RWORK,&INFO);
                
            delete[] RWORK;
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dgesvd_ ROUTINE     ////////////////////////////
        	
            dgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,
                &LDVT,WORK,&LWORK,&INFO);			
        }
        
        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////
               
        delete[] WORK;
        
        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////
        
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zgesvd_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dgesvd_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function calls a LAPACK routine to evaluate the singular
    //! values and singular vectors of a rectangular matrix.
    //!
    //! NOTE: the matrix A is overwritten with the left singular vectors
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool RightSingularValueDecomposition(
        T* A,                   //!<    Memory address of matrix for which the SVD is required
        T* VT,                  //!<    Memory address to store the right singular vectors 
        double* S,              //!<    Memory address of singular values
        const int nbrRows,      //!<    Number of rows in the matrix
        const int nbrCols)      //!<    Number of columns in the matrix
    {
        //////      SET PARAMETER VALUES FOR zgesvd_/dgesvd_ LAPACK ROUTINE     ////////
        
        char JOBU  = 'O';       	//  Setting to store the left singular vectors
                                    //  OVERWRITE the input matrix A
        char JOBVT = 'S';       	//  Setting to store the right singular vectors
                                    //  in the VT array
        int M      = nbrRows;   	//  Number of rows in the input matrix
        int N      = nbrCols;   	//  Number of columns in the input matrix
        int LDA    = std::max(1,M); //  The leading dimension of the input matrix
        T* U = 0;                   //  Not addessed with current JOBU setting
        int LDU    = M;             //  Leading dimension of the array U
        int LDVT   = N;             //  Leading dimension of the array VT 
        int LWORK  = std::max(1,2*N+2*M);
                                    //  Working memory space size
        T *WORK = new T[LWORK];     //  Working memory space
        int INFO;               	//  Diagnostic return value

        //  Change from column major to row major indexing
        //  in order to correctly pass to fortran routine	
        
        //dcmplx* B  = new dcmplx[nbrRows*nbrCols];//   Matrix to calculate the SVD for
        
        //for(int i=0;i<nbrRows;i++)
        //{
        //    for(int j=0;j<nbrCols;j++)
        //    {   
        //        B[i*nbrCols+j] = A[j*nbrRows+i];
        //    }
        //}
           
        if(utilities::is_same<T,dcmplx>::value)
        {
            double *RWORK = new double[LWORK];
                                    //  Working memory space
        
            //////      CALL LAPACK zgesvd_ ROUTINE     ////////////////////////////
        	
            zgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,
                &LDVT,WORK,&LWORK,RWORK,&INFO);
                
            delete[] RWORK;
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dgesvd_ ROUTINE     ////////////////////////////
        	
            dgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,
                &LDVT,WORK,&LWORK,&INFO);			
        }
        
        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////
               
        delete[] WORK;
        
        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////
        
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zgesvd_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dgesvd_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function calls a LAPACK routine to evaluate the singular
    //! values of a complex rectangular matrix.
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool SingularValueDecomposition(
        T* A,                   //!<    Memory address of matrix for which the SVD is required
        double* S,              //!<    Memory address of singular values
        const int nbrRows,      //!<    Number of rows in the matrix
        const int nbrCols)      //!<    Number of columns in the matrix
    {
        //////      SET PARAMETER VALUES FOR zgesvd_/dgesvd_ LAPACK ROUTINE     ////////
        
        char JOBU  = 'N';       	//   No left singular vectors are computed
        char JOBVT = 'N';       	//   No right singular vectors are computed
        int  M     = nbrRows;   	//   Number of rows in the input matrix
        int  N     = nbrCols;   	//   Number of columns in the input matrix
        //dcmplx* A  = new dcmplx[nbrRows*nbrCols];
                                    //   Matrix to calculate the SVD for
        int LDA    = nbrRows;   	//   The leading dimension of the input matrix
        T* U  = 0;                  //   Stores the left singular vectors if requested
        int LDU    = 1;             //   Leading dimension of the array U
        T* VT = 0;                  //   Stores the right singular vectors if requested
        int LDVT   = 1;             //   Leading dimension of the array VT 
        int LWORK = std::max(1,10*nbrCols);
                                    //  Working memory space size
        T *WORK = new T[LWORK];     //   Working memory space
        int INFO;               	//  Diagnostic return value
          
        //  Change from column major to row major indexing
        //  in order to correctly pass to fortran routine	
        
        //for(int i=0;i<nbrRows;i++)
        //{
        //    for(int j=0;j<nbrCols;j++)
        //    {   
        //        A[i*nbrCols+j] = matrix[j*nbrRows+i];
        //    }
        //}  
            
        if(utilities::is_same<T,dcmplx>::value)
        {
            double *RWORK = new double[LWORK];
                                    //  Working memory space
        
            //////      CALL LAPACK zgesvd_ ROUTINE     ////////////////////////////
        	
            zgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,
                &LDVT,WORK,&LWORK,RWORK,&INFO);
                
            delete[] RWORK;
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dgesvd_ ROUTINE     ////////////////////////////
        	
            dgesvd_(&JOBU,&JOBVT,&M,&N,A,&LDA,S,U,&LDU,VT,
                &LDVT,WORK,&LWORK,&INFO);			
        }	
        
        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////
               
        delete[] WORK;
        
        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////
        
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zgesvd_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dgesvd_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function calls a LAPACK routine to evaluate the eigenvalues
    //! of a real symmetric or Hermitian matrix
    //!
    //! NOTE: the input matrix will be destroyed when calling this function
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool DiagonalizeSymmetricMatrix(
        T* A,                   //!<    Memory address of matrix for which the eigenvalues are required
        double* eigenvalues, 	//!<    Memory address of eigenvalues (MUST be of length dimension)
        const int dimension,    //!<    Number of rows in the matrix
        char UPLO)              //!<    Specify wither the upper or lower triangular 
                                //!     elements are stored
    {
        //////      SET PARAMETER VALUES FOR zheevr_/dsyevr_ LAPACK ROUTINE     ////////

        char JOBZ='N';		//	Compute eigenvalues only
        char RANGE='A';		//	All eigenvalues will be found
            
        //  NOTE: Fortran indexes the columns and rows in a different way to C.
        //  We can take that into account by flipping the UPLO definition

        if('U' == UPLO)
        {
            UPLO = 'L';
        }
        else
        {
            UPLO = 'U';
        }

        int N = dimension;		//	Matrix dimension
        int LDA = N;		    //	Leading m_dimension of the matrix
        double VL;              //  Lower bound on eigenvalue interval
                                //	(Not used for RANGE='A' option)
        double VU;              //  Upper bound on eigenvalue interval
                                //	(Not used for RANGE='A' option)
        int IL;                 //  Index of the smallest eigenvalue
        int IU;                 //  Index of the largest eigenvalue
        double ABSTOL = 0;	    //	Absolute error tolerance (default set)
        int M;				    //	Total number of eigenvalues found (on output)
        T *Z = new T[N*N];      //	Memory allocation to store eigenvectors
        int LDZ = N;            //  Leading dimension of the eigenvector array
        int *ISUPPZ = new int[2*N];
         	                    //	To store eigenvalue indexes
        int LWORK = N*(30);     //  Size of LWORK allocation
        T *WORK = new T[LWORK]; //	Memory allocation for work space
        int LIWORK = N*(10);    //  Size of IWORK allocation
        int *IWORK = new int [LIWORK];	    
                                //	Memory allocation for work space
        int INFO;			    //	output flag
 
        if(utilities::is_same<T,dcmplx>::value)
        {
            //////      CALL LAPACK zheevr_ ROUTINE     ////////////////////////////

            int LRWORK = N*(30);                //  Size of RWORK allocation
            double *RWORK = new double[LRWORK];	//	Memory allocation for work space

            zheevr_(&JOBZ,&RANGE,&UPLO,&N,A,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,
                    eigenvalues,Z,&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&INFO);
            
            delete[] RWORK;
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dsyevr_ ROUTINE     ////////////////////////////

            dsyevr_(&JOBZ,&RANGE,&UPLO,&N,A,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,
                    eigenvalues,Z,&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&INFO);
        }

        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////

        delete[] WORK;
        delete[] IWORK;
        delete[] Z;
        delete[] ISUPPZ;

        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////
    
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zheevr_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dsyevr_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function calls a LAPACK routine to evaluate the eigenvalues
    //! and eigenvectors of a real symmetric or Hermitian matrix
    //!
    //! NOTE: the LAPACK routine generates memory errors if a diagonal matrix
    //! is passed to it.
    //!
    //! NOTE: the input matrix will be destroyed when calling this function
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool DiagonalizeSymmetricMatrix(
        T* A,                   //!<    A matrix of dimension N by N
        T* eigenvectors,        //!<    A memory address (N by N in size) to store
                                //!     the resulting eigenvectors
        double* eigenvalues,    //!<    The memory address (N in size) to store
                                //!     the resulting eigenvalues
        const int dimension,    //!<    The dimension N of the matrix
        const int numberToFind, //!<    Only calculate the eigenvectors associated
                                //!     with the lowest numberToFind eigenvalues
        char UPLO)              //!<    Specify wither the upper or lower triangular 
                                //!     elements are stored
    {
        //////      SET PARAMETER VALUES FOR zheevr_/dsyevr_ LAPACK ROUTINE     ////////

        char JOBZ='V';		//	Compute eigenvalues and eigenvectors
        char RANGE='I';		//	Only specified eigenvectors will be found 
            
        //  NOTE: Fortran indexes the columns and rows in a different way to C.
        //  We can take that into account by flipping the UPLO definition

        if('U' == UPLO)
        {
            UPLO = 'L';
        }
        else
        {
            UPLO = 'U';
        }

        int N = dimension;		//	Matrix dimension
        int LDA = N;		    //	Leading m_dimension of the matrix
        double VL;              //  Lower bound on eigenvalue interval
                                //	(Not used for RANGE='A' option)
        double VU;              //  Upper bound on eigenvalue interval
                                //	(Not used for RANGE='A' option)
        int IL = 1;             //  Index of the smallest eigenvalue
        int IU = numberToFind;  //  Index of the largest eigenvalue
        double ABSTOL = 0;	    //	Absolute error tolerance (default set)
        int M;				    //	Total number of eigenvalues found (on output)
        double* W = new double[N];
                                //  Working space/to contain output eigenvalues
        int LDZ = N;            //  Leading dimension of the eigenvector array
        int *ISUPPZ = new int[2*std::max(1,numberToFind)];
         	                    //	To store eigenvalue indexes
        int LWORK = N*(30);     //  Size of LWORK allocation
        T *WORK = new T[LWORK]; //	Memory allocation for work space
        int LIWORK = N*(10);    //  Size of IWORK allocation
        int *IWORK = new int[LIWORK];	    
                                //	Memory allocation for work space
        int INFO;			    //	output flag

        if(utilities::is_same<T,dcmplx>::value)
        {
            //////      CALL LAPACK zheevr_ ROUTINE     ////////////////////////////

            int LRWORK = N*(30);                //  Size of RWORK allocation
            double *RWORK = new double[LRWORK];	//	Memory allocation for work space

            zheevr_(&JOBZ,&RANGE,&UPLO,&N,A,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,
                    W,eigenvectors,&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&INFO);
            
            delete[] RWORK;
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dsyevr_ ROUTINE     ////////////////////////////

            dsyevr_(&JOBZ,&RANGE,&UPLO,&N,A,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,
                    W,eigenvectors,&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&INFO);
        }

        for(int i=0;i<numberToFind;i++)
        {
            eigenvalues[i] = W[i];
        }

        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////

        delete[] WORK;
        delete[] IWORK;
        delete[] ISUPPZ;
        delete[] W;

        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////
    
        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zheevr_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dsyevr_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief	Given a Hermitian/symmetric input matrix M, this function calculates
    //! M^power where ``power'' is a real values exponent
    //!
    //! The routine works by diagonalizing the matrix as U^dagger D U. Then the 
    //! matrixpower is defined by:
    //!
    //! M^power = U^dagger D^power U
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    void SymmetricMatrixPower(
        T* A,                    //!<   Matrix to be operated on
        const int dimension,     //!<   Matrix dimension
        const double power)      //!<   Power to which the matrix is to be raised
    {
        //  Allocate memory space to store all the eigenvectors
        
        T* eigenvectors = new T[dimension*dimension];
        double eigenvalues[dimension];

        //  Diagonalize the matrix and return all eigenvectors;
        
        //  Make a copy of the matrix to diagonalize
        
        T *A1 = new T[dimension*dimension];
        
        memcpy(A1,A,dimension*dimension*sizeof(T));
        
        DiagonalizeSymmetricMatrix<T>(A1,eigenvectors,eigenvalues,dimension,dimension,'U');
        
        delete[] A1;
        
        //std::cout<<"\n\tEIGENVALUES: "<<std::endl;
        //for(int i=0;i<dimension;i++)
        //{
        //    std::cout<<eigenvalues[i]<<std::endl;
        //}
        
        //std::cout<<"\n\tEIGENVECTOR MATRIX: "<<std::endl;
        
        //for(int i=0;i<dimension;i++)
        //{
        //    for(int j=0;j<dimension;j++)
        //    {
        //        std::cout<<eigenvectors[i*dimension+j]<<"\t";
        //   }    
        //    std::cout<<std::endl;
        //}
        
        //  Take the regular exponent of the eigenvalues
        
        for(int i=0;i<dimension;i++)
        {
            eigenvalues[i] = pow(eigenvalues[i],power);
        }

        //  Reconstruct the matrix from U^dagger D^power U
        //  In the process, overwrite the original matrix

        T *p_m = A;
        
        for(int i=0;i<dimension;i++)
        {
            for(int j=0;j<dimension;j++,p_m++)
            {
                *(p_m) = 0;
                
                double *p_D = eigenvalues;
                T *p_U = eigenvectors + j;    
                T *p_Udagger = eigenvectors + i;

                for(int k=0;k<dimension;k++,p_D++,p_Udagger+=dimension,p_U+=dimension)
                { 
                    //std::cout<<"\n\tBEGIN: i = "<<i<<"j = "<<j<<" k = "<<k<<std::endl;
                
                    //  Implement H_ij = sum_k U^dagger_ik D_kk U_kj 
                    //  so H_ij = sum_k U*_ki D_kk U_kj
            
                    if(utilities::is_same<T,dcmplx>::value)
                    {
                        *(p_m) += std::conj(*(p_Udagger))**(p_D)**(p_U);
                    }
                    else if(utilities::is_same<T,double>::value)
                    {
                        *(p_m) += *(p_Udagger) * *(p_D) * *(p_U);
                    }
                    
                    //std::cout<<"\n\tEND: i = "<<i<<"j = "<<j<<" k = "<<k<<std::endl;
                }
            }
     
        }

        delete[] eigenvectors;

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function calls a LAPACK diagonalization routine to diagonalize
    //! a non-symmetric real matrix.
    //!
    //! NOTE: a transposed copy of the input matrix is taken (so the input matrix 
    //! is not destroyed by this function)
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool DiagonalizeGeneralMatrix(
        T* A,                       //!<    Pointer to memory address of the matrix
        T* leftVectors,             //!<    Pointer to where the eigenvectors are stored
        T* rightVectors,            //!<    Pointer to where the eigenvectors are stored
        dcmplx* eigenvalues,        //!<    Pointer to where the eigenvalues are stored
        const int dimension)        //!<    Dimension of the matrix
    {
        //////      SET PARAMETER VALUES FOR dgeev_/zgeev_ LAPACK ROUTINE     ////////
    
        char JOBVL = 'V';       //  Option to evaluate left eigenvectors
        char JOBVR = 'V';       //  Option to evaluate right eigenvectors
        int N = dimension;      //  Dimension of the matrix
        int LDA = N;            //  Leading dimension of the matrix        
        int LDVL = N;           //  Dimension of VL array
        int LDVR = N;           //  Dimension of VR array
        int LWORK = 20*N;       //  Dimension of working memory array
        T *WORK = new T[LWORK]; //  Working memory array
        int INFO;               //  Exception flag
        
        //  Change from column major to row major indexing
        //  in order to correctly pass to fortran routine	
        
        //dcmplx* B  = new dcmplx[nbrRows*nbrCols];//   Matrix to calculate the SVD for
        
        //for(int i=0;i<nbrRows;i++)
        //{
        //    for(int j=0;j<nbrCols;j++)
        //    {   
        //        B[i*nbrCols+j] = A[j*nbrRows+i];
        //    }
        //}
        
        if(utilities::is_same<T,dcmplx>::value)
        {
            zgeev_(&JOBVL,&JOBVR,&N,A,&LDA,eigenvalues,leftVectors,&LDVL,rightVectors,
                   &LDVR,WORK,&LWORK,&INFO);
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dgeev_ ROUTINE     ////////////////////////////
            
            double *WR = new double[N]; //  To contain real parts of computed eigenvalues
            double *WI = new double[N]; //  To contain imaginary parts of computed eigenvalues
            
            dgeev_(&JOBVL,&JOBVR,&N,A,&LDA,WR,WI,leftVectors,&LDVL,rightVectors,
                   &LDVR,WORK,&LWORK,&INFO);
                   
            for(int k=0;k<dimension;k++)
            {
                eigenvalues[k] = dcmplx(WR[k],WI[k]);
                
                //std::cout<<"\tEIGENVALUES: "<<eigenvalues[k]<<std::endl;
            }
            
            delete[] WR;
            delete[] WI;
        }
        
        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////

        delete[] WORK;

        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////

        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zgeev_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dgeev_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function calls a LAPACK diagonalization routine to diagonalize
    //! a non-symmetric real matrix.
    //!
    //! NOTE: a transposed copy of the input matrix is taken (so the input matrix 
    //! is not destroyed by this function)
    //!
    //! \return true if there was an error, false otherwise
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    bool DiagonalizeGeneralMatrix(
        T* A,                       //!<    Pointer to memory address of the matrix
        dcmplx* eigenvalues,        //!<    Pointer to where the eigenvalues are stored
        const int dimension)        //!<    Dimension of the matrix
    {
        //////      SET PARAMETER VALUES FOR dgeev_ LAPACK ROUTINE     ////////
    
        char JOBVL = 'N';       //  Do not evaluate left eigenvectors
        char JOBVR = 'N';       //  Do not evaluate right eigenvectors
        int N = dimension;      //  Dimension of the matrix
        int LDA = N;            //  Leading dimension of the matrix        
        int LDVL = N;           //  Dimension of VL array
        T* VL = 0;              //  Not addressed
        int LDVR = N;           //  Dimension of VR array
        T* VR = 0;              //  Not addressed
        int LWORK = 20*N;       //  Dimension of working memory array
        T *WORK = new T[LWORK]; //  Working memory array
        int INFO;               //  Exception flag
        
        //  Change from column major to row major indexing
        //  in order to correctly pass to fortran routine	
        
        //dcmplx* B  = new dcmplx[nbrRows*nbrCols];//   Matrix to calculate the SVD for
        
        //for(int i=0;i<nbrRows;i++)
        //{
        //    for(int j=0;j<nbrCols;j++)
        //    {   
        //        B[i*nbrCols+j] = A[j*nbrRows+i];
        //    }
        //}
        
        if(utilities::is_same<T,dcmplx>::value)
        {
            zgeev_(&JOBVL,&JOBVR,&N,A,&LDA,eigenvalues,VL,&LDVL,VR,
                   &LDVR,WORK,&LWORK,&INFO);
        }
        else if(utilities::is_same<T,double>::value)
        {
            //////      CALL LAPACK dgeev_ ROUTINE     ////////////////////////////
            
            double *WR = new double[N]; //  To contain real parts of computed eigenvalues
            double *WI = new double[N]; //  To contain imaginary parts of computed eigenvalues
            
            dgeev_(&JOBVL,&JOBVR,&N,A,&LDA,WR,WI,VL,&LDVL,VR,
                   &LDVR,WORK,&LWORK,&INFO);
                   
            for(int k=0;k<dimension;k++)
            {
                eigenvalues[k] = dcmplx(WR[k],WI[k]);
                
                //std::cout<<"\tEIGENVALUES: "<<eigenvalues[k]<<std::endl;
            }
            
            delete[] WR;
            delete[] WI;
        }
        
        //////      CLEAR UP MEMORY ALLOCATION      ////////////////////////////

        delete[] WORK;

        //////      CHECK FOR ERROR MESSAGES        ////////////////////////////

        if(INFO!=0)
        {
            if(utilities::is_same<T,dcmplx>::value)
            {
                std::cerr<<"\n\tERROR with zgeev_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            else if(utilities::is_same<T,double>::value)
            {
                std::cerr<<"\n\tERROR with dgeev_ algorithm : INFO returned "<<INFO<<std::endl;
            }
            
            return true;
        }

        return false;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace linearAlgebra

}   //  End namespace utilities

#endif

