////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file
//!		This is my header file for sparse linear algebra functions. 
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

#ifndef _SPARSE_LINEAR_ALGEBRA_HPP_INCLUDED_
#define _SPARSE_LINEAR_ALGEBRA_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../general/dcmplx_type_def.hpp" 
#include "../data_structures/crs_sparse_matrix.hpp" 
#include "../general/cout_tools.hpp" 
#include "../wrappers/mpi_wrapper.hpp"
#include "../general/template_tools.hpp"
#include <vector>
#if _DEBUG_
#include "../general/debug.hpp"
#endif

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_INTEL_MKL_
//////      IMPORT INTEL MKL (P)BLAS ROUTINES FOR   ////////////////////////////
//////      SPARSE MATRIX-VECTOR MULTIPLICATION     ////////////////////////////
extern "C"
{
    void mkl_zcsrmv(char* TRANSA, int* M, int* K, dcmplx* ALPHA, char* MATDESCRA,
                    dcmplx* VAL, int* INDX, int* PNTRB, int* PNTRE, dcmplx* X, 
                    dcmplx* BETA, dcmplx* Y, int* LEN_MATDESCRA);           

    void mkl_dcsrmv(char* TRANSA, int* M, int* K, double* ALPHA, char* MATDESCRA,
                    double* VAL, int* INDX, int* PNTRB, int* PNTRE, double* X, 
                    double* BETA, double* Y, int* LEN_MATDESCRA);           
}
//////      Dummy declarations to allow for different       ////////////////////
//////      template arguments. Functions not referenced.   ////////////////////
void mkl_zcsrmv(char* TRANSA, int* M, int* K, double* ALPHA, char* MATDESCRA,
                double* VAL, int* INDX, int* PNTRB, int* PNTRE, double* X, 
                double* BETA, double* Y, int* LEN_MATDESCRA);
                    
void mkl_dcsrmv(char* TRANSA, int* M, int* K, dcmplx* ALPHA, char* MATDESCRA,
                dcmplx* VAL, int* INDX, int* PNTRB, int* PNTRE, dcmplx* X, 
                dcmplx* BETA, dcmplx* Y, int* LEN_MATDESCRA);                    
#endif
////////////////////////////////////////////////////////////////////////////////

namespace utilities
{
namespace linearAlgebra
{
    //////////////////////////////////////////////////////////////////////////////////   
    //! \brief A function to perform matrix-vector multiplication on sparse
    //! type matrices and vectors: y = alpha*A*x+beta*x where A is a Symmetric
    //! or Hermitian sparse matrix
    //!
    //! This function calls an MKL BLAS routine to perform Hermitian 
    //! matrix-vector  multiplication. If MKL is not available, use a naive 
    //! and slow back up method.
    //////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    void SymmetricMatrixVectorMultiply(
        T* x,           //!<    Input vector
        T* y,           //!<    Output vector
        int dim,        //!<    Dimension of the vectors/square matrix
        utilities::CrsSparseMatrix<T>* matrix,  
                        //!<	Address of the sparse CRS format matrix
        char* UPLO,     //!<    Specify 'U' if upper triangular is stored, 
                        //!<    of 'L' if lower triangular is stored
        T ALPHA,        //!<    Matrix scale factor
        T BETA)         //!<    Additional vector scale factor
    {
        ////////////////////////////////////////////////////////////////////////////////
        #if _ENABLE_INTEL_MKL_
        //////      Declare variables to pass to mkl blas routine
        //  void mkl_?csrmv(char* TRANSA,int* M,int* K,T* ALPHA,char* MATDESCRA,
        //    T* VAL,int* INDX,int* PNTRB,int* PNTRE,T* X, 
        //    T* BETA,T* Y,int* LEN_MATDESCRA)
        char TRANSA;        //  if TRANSA = 'N' then y = alpha*A*x+beta*y
                            //  if TRANSA = 'T' then y = alpha*A'*x+beta*y
        int M;              //  Number of rows in the matrix A
        int K;              //  Number of columns in the matrix A
        char MATDESCRA[6];  //  Matrix description (see Intel's website 
                            //  for possible values)
        T* VAL;             //  Array containing non-zero elements                    
        int* INDX;          //  Array containing column indices
        int* PNTRB;         //  Array of dimension M such that PNTRB[i] - PNTRB[0] 
                            //  is the first index of row i
        int* PNTRE;         //  Array of dimension M such that PNTRE[i] - PNTRB[0] - 1
                            //  Is the last index of row i
        T* X;               //  Input vector             
        T* Y;               //  Output vector  
        int LEN_MATDESCRA;  //  Additional argument to pass the length of the
                            //  MATDESCRA character array                      
        //////      Set parameter values
        TRANSA = 'N';           //  No transpose
        M = dim;
        K = dim;
        MATDESCRA[0] =  'H';    //  Hermitian matrix
        MATDESCRA[1] =  *UPLO;  //  Upper or lower triangular specified          
        MATDESCRA[2] =  'N';    //  Not a unit matrix
        MATDESCRA[3] =  'C';    //  Use C-syle matrix indexing
        VAL  = matrix->GetData();
        INDX = matrix->GetCols();
        PNTRB = matrix->GetRowStarts();
        PNTRE = matrix->GetRowStarts()+1;
        X = x;
        Y = y;
        LEN_MATDESCRA = 6;
        //  Copy the x into the y vector. This allows us to calculate
        //  y = alpha Ax + beta*x (since the mkl function actually does 
        //  y = alpha Ax + beta*y)
        if(BETA!=0.0)
        {
            memcpy(Y, X, (long int)dim*sizeof(T));
        }
        if(utilities::is_same<T, dcmplx>::value)
        {
            mkl_zcsrmv(&TRANSA, &M, &K, &ALPHA, &MATDESCRA[0], VAL, INDX, PNTRB, PNTRE, X, &BETA, Y, &LEN_MATDESCRA);
        }
        else if(utilities::is_same<T, double>::value)
        {
            mkl_dcsrmv(&TRANSA, &M, &K, &ALPHA, &MATDESCRA[0], VAL, INDX, PNTRB, PNTRE, X, &BETA, Y, &LEN_MATDESCRA);
        }
        #else
        //  Naive (and slow) backup implementation
        for(long int i=0; i<dim; ++i)
        {
            y[i] = BETA*x[i]; 
        }
        int* p_rowStarts = matrix->GetRowStarts();
        int* p_cols      = matrix->GetCols();
        T* p_data        = matrix->GetData();
        for(int i=0; i<matrix->GetNbrRows(); ++i)
        {
            for(int j = p_rowStarts[i]; j<p_rowStarts[i+1]; ++j)
            {
                int tmpCol   = p_cols[j];
                T tmpElement = p_data[j];
                y[i] += ALPHA*tmpElement*x[tmpCol];
                //  Take into account that our matrix is Hermitian
                if(i !=tmpCol)
                {
                    y[tmpCol] += ALPHA*std::conj(tmpElement)*x[i];
                }
            }
        }
        #endif
        ////////////////////////////////////////////////////////////////////////////////
        return;
    }

    //////////////////////////////////////////////////////////////////////////////////   
    //! \brief A function to perform a parallel matrix-vector multiplication on sparse
    //! type matrices and vectors: y = A*x where A is a symmetric or
    //! Hermitian sparse matrix
    //!
    //! This function is optimized to try to reduce inter-node communication.
    //!
    //! NOTE 1: The function expects that the input and output vectors will be
    //! distributed over a given number of nodes. The dataDistribution
    //! argument specifies the dimension of the part of the input/output
    //! vector stored on the each node (and is identical on all nodes). 
    //!
    //! NOTE 2: In order for this function to work properly, the sparse matrix
    //! information must be updated with the PrepareForParallelMultiplication
    //! function.
    //!
    //! The picture to have in mind is as follows. We solve Y = A X, where
    //! the memory distribution of the vectors/matrix is as:
    //!
    //!         +------------+
    //! Y0      | A0         |     X0   }   dataDistribution[0]
    //!         |   +--------+
    //! Y1      |   | A1     |     X1   }   dataDistribution[1]
    //!         |   |   +----+
    //! Y2      |   |   | A2 |     X2   }   dataDistribution [2]
    //!
    //! Where the first dataDistribution[0] rows of the matrix are stored on 
    //! node 0, and so on. Note that the Hermitian structure of the matrix 
    //! means that the first n rows also give the first n columns.
    //! 
    //! The matrix A0 can be thought of as a NxN matrix (N the total dimension)
    //! with zeros where A1, A2 are labelled. A1 is a (N-dataDistribution[0])x 
    //! (N-dataDistribution[0]) matrix with zeros where A2 is labelled, and so on.
    //! Appropriate restructuring is done by PrepareForParallelMultiplication.
    //!
    //! The matrix-vector multiplication can be done like this:
    //!
    //! Y = A0.[X0,X1,X2] + A1.[X1,X2] + A2.[X2]
    //!
    //! Where e.g. [X0,X1,X2] is a vector combining X0,X1,X2 and where the 
    //! output vectors of each matrix multiplication are appropriately padded 
    //! with zeros when summed. 
    //! 
    //! Assuming that the As take up a lot of memory, we focus on passing
    //! the minimal amount of X around. Note also that the result in e.g Y0
    //! does not depend on A1.[X1,X2] or A2.[X2] etc, so we can also 
    //! minimise the amount of Y passed around.
    ////////////////////////////////////////////////////////////////////////////////// 
    template <typename T>
    void ParallelSymmetricMatrixVectorMultiply(
        T* x,                               //!<    Input vector on the current node
        T* y,                               //!<    Output vector on the current node
        std::vector<int>* dataDistribution, //!<    A vector containing a list
                                            //!     of the number of matrix rows
                                            //!     stored on each node
        utilities::CrsSparseMatrix<T>* matrix,  
                                            //!<    Address of the sparse CRS format matrix
        char* UPLO,                         //!<    Specify 'U' if upper triangular is stored,
                                            //!<    of 'L' if lower triangular is stored
        T ALPHA,                            //!<    Matrix scale factor
        T BETA,                             //!<    Additional vector scale factor
        std::vector<MPI_Comm>* commGroups,  //!<    List of processor groups
                                            //!     required for minimal
                                            //!     vector communication
        T* inputVectorBuffer,               //!<    Buffer for MPI distribution of X vector
        T* outputVectorBuffer,              //!<    Buffer for MPI accumulation of Y vector
        const int bufferDim,                //!<    Size of the input/output buffers
        const MpiWrapper& mpi)              //!<    Instance of the MPI wrapper class
    {
        //  The vector scale factor BETA is only used on the master node, 
        //  and should be set to zero everywhere else
        if(0 != mpi.m_id)
        {
            BETA = 0.0;
        }
        //  To implement this, we define a set of M-1 processor groups, 
        //  where M is the number of processors. The 1st group contains
        //  all M nodes, the second group contains M-1 nodes, etc. 
        //  The first step is to put exactly the required amount of the input
        //  vector on each node. Node i requires only Xi, Xi+1 and so on,
        //  but not Xi-1 etc.
        int cumulativeBcast = 0;      
        for(int i = 0; i<mpi.m_nbrProcs; ++i)
        { 
            if(mpi.m_id == i)
            {
                memcpy(inputVectorBuffer, x, (*dataDistribution)[i]*sizeof(T));
            }
            if(0 != i)
            {
                //  Broadcast the contents of the buffer to the preceding i nodes 
                //  only, by splitting up the MPI communicator into 2 groups -
                //  one for the preceding nodes, one for the rest
                int whichGroup = (i>=mpi.m_id ? 0 : 1);
                //  Only call the Bcast function on one of the processor
                //  groups - this saves a factor of 2 in communication time!
                if(0 == whichGroup)
                {
                    MPI_Bcast(inputVectorBuffer+cumulativeBcast, (*dataDistribution)[i], mpi.GetType<T>(), i, (*commGroups)[i]);
                }
            }
            if(mpi.m_id <= i)
            {
                cumulativeBcast += (*dataDistribution)[i];
            }
        }
        //  Now inputVectorBuffer on node i stores [Xi,Xi+1,...]
        //  Next we can perform a serial matrix-vector multiplication
        //  i.e. the form Ai.[Xi,Xi+1,...]. The result will be stored in
        //  the outputVectorBuffer.
        SymmetricMatrixVectorMultiply(inputVectorBuffer, outputVectorBuffer, bufferDim, matrix, UPLO, ALPHA, BETA);
        //  Finally we re-distribute the outputVectorBuffer to generate 
        //  the correct Y vector on each node. This step performs the
        //  sum e.g. Y = A1.[X1,X2,X3] + A2.[X2,X3] + A3.[X3]
        //  Using an MPI reduce algorithm
        int cumulativeReduce = 0;
        for(int i = 0; i<mpi.m_nbrProcs; ++i)
        { 
            if(0 != i)
            {
                int whichGroup = (i>=mpi.m_id ? 0 : 1);
                //  Only call the Reduce function on one of the processor
                //  groups - this saving a factor of 2 in communication time!
                if(0 == whichGroup)
                {
                    //  Note that the inputVectorBuffer is recycled 
                    //  as working memory for the mpi.Reduce function
                    MPI_Status status;
                    mpi.Reduce<T>(outputVectorBuffer+cumulativeReduce,inputVectorBuffer,
                                  (*dataDistribution)[i], i, (*commGroups)[i], status);
                }
            }
            if(mpi.m_id == i)
            {
                memcpy(y, outputVectorBuffer, (*dataDistribution)[i]*sizeof(T));
            }
            if(mpi.m_id <= i)
            {
                cumulativeReduce += (*dataDistribution)[i];
            }
        }
        return;
    }
}   //  End namespace linearAlgebra 
}   //  End namespace utilities

#endif
