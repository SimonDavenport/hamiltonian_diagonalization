////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file
//!     This file provides wrappers for several ARPACK and PARPACK routines.
//!     See here http://www.caam.rice.edu/software/ARPACK/ for documentation.
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

#ifndef _ARPACK_WRAPPER_HPP_INCLUDED_
#define _ARPACK_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../general/dcmplx_type_def.hpp"
#include "../general/cout_tools.hpp"
#include "../wrappers/mpi_wrapper.hpp"
#include <vector>
#include <functional> 
#if _DEBUG_
#include "../general/debug.hpp"
#endif

///////		Import ARPACK routines to diagonalize       ////////////////////////
//////      sparse Hermitian matrices                   ////////////////////////
extern "C"
{
////////////////////////////////////////////////////////////////////////////////
//! \brief znaupd_ is an ARPACK routine for diagonalizing Hermitian matrices.
//! znaupd_ is a parallel architecture version
//!
//! "Reverse communication interface for the Implicitly Restarted Arnoldi
//! iteration. This is intended to be used to find a few eigenpairs of a 
//! complex linear operator OP with respect to a semi-inner product defined 
//! by a hermitian positive semi-definite real matrix B. B may be the identity 
//! matrix."
//!
//! http://www.mathkeisan.com/usersguide/man/znaupd.html
//!
//! http://www.mathkeisan.com/usersguide/man/pznaupd.html
////////////////////////////////////////////////////////////////////////////////
void znaupd_(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
             dcmplx* RESID, int* NCV, dcmplx* V, int* LDV, int* IPARAM, int* IPNTR,
             dcmplx* WORKD, dcmplx* WORKL, int* LWORKL, double* RWORK, int* INFO, 
             int* LEN_BMAT, int* LEN_WHICH);
            
void pznaupd_(int* COMM, int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
              dcmplx* RESID, int* NCV, dcmplx* V, int* LDV, int* IPARAM, int* IPNTR,
              dcmplx* WORKD, dcmplx* WORKL, int* LWORKL, double* RWORK, int* INFO, 
              int* LEN_BMAT, int* LEN_WHICH);
////////////////////////////////////////////////////////////////////////////////
//! \brief dnaupd_ is an ARPACK routine for diagonalizing Symmetric matrices.
//! pdnaupd_ is a parallel architecture version
//!
//! http://www.mathkeisan.com/usersguide/man/dnaupd.html
//!
//! http://www.mathkeisan.com/usersguide/man/pdnaupd.html
////////////////////////////////////////////////////////////////////////////////
void dnaupd_(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
             double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR,
             double* WORKD, double* WORKL, int* LWORKL, int* INFO, int* LEN_BMAT, int* LEN_WHICH);

void pdnaupd_(int* COMM, int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
              double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR,
              double* WORKD, double* WORKL, int* LWORKL, int* INFO, int* LEN_BMAT, int* LEN_WHICH);
////////////////////////////////////////////////////////////////////////////////
//! \brief zneupd_ is an ARPACK routine for diagonalizing Hermitian matrices
//! zneupd_ is a parallel architecture version
//!
//! "Postprocessing routine for large-scale complex eigenvalue calculation."
//!
//! NOTE: ERROR IN DOCUMENATION ON ARGUMENT ORDER
//!       
//! http://www.mathkeisan.com/UsersGuide/man/zneupd.html
//!
//! http://www.mathkeisan.com/usersguide/man/pzneupd.html
////////////////////////////////////////////////////////////////////////////////          
void zneupd_(int* RVEC, char* HOWMNY, int* SELECT, dcmplx* D, dcmplx* Z, int* LDZ,
             dcmplx* SIGMA, dcmplx* WORKEV, char* BMAT, int* N, char* WHICH, int* NEV,
             double* TOL, dcmplx* RESID, int* NCV, dcmplx* V, int* LDV, int* IPARAM,
             int* IPNTR, dcmplx* WORKD, dcmplx* WORKL, int* LWORKL, double* RWORK,
             int* INFO, int* LEN_HOWMNY, int* LEN_BMAT, int* LEN_WHICH);

void pzneupd_(int* COMM, int* RVEC, char* HOWMNY, int* SELECT, dcmplx* D, dcmplx* Z, 
              int* LDZ, dcmplx* SIGMA, dcmplx* WORKEV, char* BMAT, int* N, char* WHICH,
              int* NEV, double* TOL, dcmplx* RESID, int* NCV, dcmplx* V, int* LDV, 
              int* IPARAM, int* IPNTR, dcmplx* WORKD, dcmplx* WORKL, int* LWORKL, 
              double* RWORK, int* INFO, int* LEN_HOWMNY, int* LEN_BMAT, int* LEN_WHICH);
////////////////////////////////////////////////////////////////////////////////
//! \brief dneupd_ is an ARPACK routine for diagonalizing Hermitian matrices
//! dneupd_ is a parallel architecture version
//!
//! NOTE: ERROR IN DOCUMENATION ON ARGUMENT ORDER
//!       
//! http://www.mathkeisan.com/UsersGuide/man/dneupd.html
//!
//! http://www.mathkeisan.com/usersguide/man/pdneupd.html
////////////////////////////////////////////////////////////////////////////////
void dneupd_(int* RVEC, char* HOWMNY, int* SELECT, double* DR, double* DI, double* Z,
             int* LDZ, double* SIGMAR, double* SIGMAI, double* WORKEV, char* BMAT, int* N,
             char* WHICH, int* NEV, double* TOL, double* RESID, int* NCV, double* V,
             int* LDV, int* IPARAM, int* IPNTR, double* WORKD, double* WORKL, int* LWORKL,
             int* INFO, int* LEN_HOWMNY, int* LEN_BMAT, int* LEN_WHICH);

void pdneupd_(int* COMM, int* RVEC, char* HOWMNY, int* SELECT, double* DR, double* DI,
             double* Z, int* LDZ, double* SIGMAR, double* SIGMAI, double* WORKEV,
             char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, double* RESID,
             int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR, double* WORKD,
             double* WORKL, int* LWORKL, int* INFO, int* LEN_HOWMNY, int* LEN_BMAT,
             int* LEN_WHICH);                
}
//////      Dummy declarations to allow for different       ////////////////////
//////      template arguments. Functions not referenced.   ////////////////////
void znaupd_(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
             double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR,
             double* WORKD, double* WORKL, int* LWORKL, double* RWORK, int* INFO, 
             int* LEN_BMAT, int* LEN_WHICH);
            
void pznaupd_(int* COMM, int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
              double* RESID, int* NCV, double* V, int* LDV, int* IPARAM, int* IPNTR,
              double* WORKD, double* WORKL, int* LWORKL, double* RWORK, int* INFO, 
              int* LEN_BMAT, int* LEN_WHICH);

void dnaupd_(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
             dcmplx* RESID, int* NCV, dcmplx* V, int* LDV, int* IPARAM, int* IPNTR,
             dcmplx* WORKD, dcmplx* WORKL, int* LWORKL, int* INFO, int* LEN_BMAT, int* LEN_WHICH);
            
void pdnaupd_(int* COMM, int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, double* TOL,
              dcmplx* RESID, int* NCV, dcmplx* V, int* LDV, int* IPARAM, int* IPNTR,
              dcmplx* WORKD, dcmplx* WORKL, int* LWORKL, int* INFO, int* LEN_BMAT, int* LEN_WHICH);
           
void zneupd_(int* RVEC, char* HOWMNY, int* SELECT, double* D, double* Z, int* LDZ,
             double* SIGMA, double* WORKEV, char* BMAT, int* N, char* WHICH, int* NEV,
             double* TOL, double* RESID, int* NCV, double* V, int* LDV, int* IPARAM,
             int* IPNTR, double* WORKD, double* WORKL, int* LWORKL, double* RWORK,
             int* INFO, int* LEN_HOWMNY, int* LEN_BMAT, int* LEN_WHICH);

void pzneupd_(int* COMM, int* RVEC, char* HOWMNY, int* SELECT, double* D, double* Z, int* LDZ,
              double* SIGMA, double* WORKEV, char* BMAT, int* N, char* WHICH, int* NEV,
              double* TOL, double* RESID, int* NCV, double* V, int* LDV, int* IPARAM,
              int* IPNTR, double* WORKD, double* WORKL, int* LWORKL, double* RWORK,
              int* INFO, int* LEN_HOWMNY, int* LEN_BMAT, int* LEN_WHICH);     

void dneupd_(int* RVEC, char* HOWMNY, int* SELECT, dcmplx* DR, dcmplx* DI, dcmplx* Z,
             int* LDZ, dcmplx* SIGMAR, dcmplx* SIGMAI, dcmplx* WORKEV, char* BMAT, int* N,
             char* WHICH, int* NEV, double* TOL, dcmplx* RESID, int* NCV, dcmplx* V,
             int* LDV, int* IPARAM, int* IPNTR, dcmplx* WORKD, dcmplx* WORKL, int* LWORKL,
             int* INFO, int* LEN_HOWMNY, int* LEN_BMAT, int* LEN_WHICH);

void pdneupd_(int* COMM, int* RVEC, char* HOWMNY, int* SELECT, dcmplx* DR, dcmplx* DI,
              dcmplx* Z, int* LDZ, dcmplx* SIGMAR, dcmplx* SIGMAI, dcmplx* WORKEV,
              char* BMAT, int* N, char* WHICH, int* NEV, double* TOL, dcmplx* RESID,
              int* NCV, dcmplx* V, int* LDV, int* IPARAM, int* IPNTR, dcmplx* WORKD,
              dcmplx* WORKL, int* LWORKL, int* INFO, int* LEN_HOWMNY, int* LEN_BMAT,
              int* LEN_WHICH);  

namespace utilities
{
    namespace linearAlgebra
    {
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A wrapper class for ARPACK's znaupd_/dnaupd_ and zneupd_/dneupd_
    //! functions (or parallel versions in PARPACK). 
    //!
    //! This class provides an interface to set the parameters up as required
    //! and call the ARPACK routines without having to specify a long list of 
    //! arguments or complicated sets of function calls. 
    //!
    //! See Arpack users guide for more information.
    ////////////////////////////////////////////////////////////////////////////////
    
    class ArpackWrapper
    {
        private:
        static const uint64_t memoryQueryLimit = 0x100000000;     
        //!<    Set 4GB max memory to store any given vector/matrix
        char m_typeOfProblem;           //!<   'I' sets a standard eigenvalue problem, 
                                        //!     'G' sets a generalized eigenvalue problem
        char m_eigenvalueMagnitude;     //!<    'S' sets smallest, 'L' sets largest
                                        //!     in set of eigenvalues to be found
        char m_whatMagnitude;           //!<    m_eigenvalueMagnitude refers to 'R'
                                        //!     for real part, 'I' for imag part or
                                        //!     'M' for absolute value
        char m_upperOrLowerMatrix;      //!<    Flag to state what part of the 
                                        //!     Hermitian matrix is stored
        //////      PROVIDE AN INTERFACE FOR THE PARAMETERS HIDDEN IN IPARAM        /////   
        int m_maxIterations;            //!<    Max number of iterations (IPARAM[2])
        static const int m_blockSize=1; //!<    Recurrence block size (IPARAM[3])
                                        //!     Only currently allowed to be 1
        //////      PROVIDE INTERFACT TO SET THE STARTING ITERATION VECTORS        /////
        int m_provideInitial;           //!<    Flag set to > 0 if the initial vector is 
                                        //!     provided - otherwise if set to 0, the 
                                        //!     initial vector is chosen at random
        bool m_storeFinalVector;        //!<    Store the final iterated vector
                                        //!     in a file, in order to pass to the 
                                        //!     initial vector of the next run
        std::string m_initialVectorFile;//!<    File where the initial vector
                                        //!     is stored/retrieved
        std::string m_finalVectorFile;  //!<    File where the final vector
                                        //!     is stored/retrieved        
        //////      DEFINE CLASS CONSTRUCTOR AND OTHER FUNCTIONS    ///////////////////
        public:
        ArpackWrapper();
        ArpackWrapper(const ArpackWrapper& other);
        ArpackWrapper& operator=(const ArpackWrapper& other);
        void MpiSync(const int nodeId, const MpiWrapper& mpi);
        //////      Public interface to edit user defined parameters    ////////
        void SetProblemType(const char typeOfProblem);
        void SetEigenvalueMagnitude(const char eigenvalueMagnitude);
        void SetWhatMagnitude(const char whatMagnitude);
        void SetUpperOrLower(const char upperOrLower);
        char GetUpperOrLower() const;
        void SetMaxIterations(const int maxIterations);
        void UseInitialVector();
        void SetInitialVectorFile(const std::string initialVectorFile);
        void StoreFinalVector();
        void SetFinalVectorFile(const std::string finalVectorFile);

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief This function calls an ARPACK diagonalization routine to obtain a
        //!	approximate eigenvalues/eigenvectors of a complex hermitian matrix
        //! or real symmetric matrix.
        //!
        //! Returns selected Ritz values/vectors (approximate eigenvalues/eigenvectors) 
        //! to the memory address specified by the eigenvalues/eigenvectors variables. 
        //!
        //! NOTE: on shift and invert mode
        //!
        //! Instead of solving Ax = lambda x we solve (A-shift*I)^-1 = nu x
        //! Then use that nu = 1/(lambda - shift)
        //!
        //! This approach reduces the number of ARPACK iterations required to 
        //! converge to low lying eigenvalues of the original spectrum, at the 
        //! expense of a more computationally demanding linear solve step at
        //! each iteration and an initial Cholesky decomposition of the 
        //! matrix (which can take a lot of memory if the decomposition is at 
        //! lot less sparse than the original matrix.) The shift value must be   
        //! chosen close to the lowest lying eigenvalue for the method to work.
        ////////////////////////////////////////////////////////////////////////////////
        template <typename T>
        bool DiagonalizeSymmetricMatrix(
        const std::function<void(T* inVector, T* outVector)>& MatrixVectorMultiply,
                                    //!<    Function to perform matrix-vector operation
        T* eigenvectors,            //!<    A memory address (N by M in size) to store
                                    //!     the resulting eigenvectors
        double* eigenvalues,        //!<    The memory address (N in size) to store
                                    //!     the resulting eigenvalues
        const int nodeDim,          //!<    Dimension of the matrix/eigenvectors 
                                    //!     stored on each node
        const int fullDim,          //!<    Full dimension of the matrix
        const int nbrEigenvalues)   //!<    Number of eigenvalues requested
        {
            //////      A DESCRIPTION OF ALL ?naupd_ ARGUMENTS     ///////////////////

            int IDO;            //!<  Must be set to zero initially,
                                //!   afterwards specifies the type of operation performed
                                //!   and value 99 tells us the iterative procedure has finished
            char BMAT;          //!<  'I' sets a standard eigenvalue problem, 'G' sets a 
                                //!  generalized eigenvalue problem instead
            int N;              //!<  Dimension of the eigenvalue problem
            char WHICH[2];      //!<  Option for range of eigenvalue returned:
                                //!  'LM' -> eigenvalues of largest magnitude
                                //!  'SM' -> eigenvalues of smallest magnitude
                                //!  'LR' -> eigenvalues of largest real part
                                //!  'SR' -> eigenvalues of smallest real part
                                //!  'LI' -> eigenvalues of largest imaginary part
                                //!  'SI' -> eigenvalues of smallest imaginary part
            int NEV;            //!<   Number of eigenvalues/eigenvectors to find
            double TOL;         //!<  Stopping tolerance (defaults to machine precision)
            T* RESID;           //!<  A dimension-N vector containing the residuals
            int NCV;            //!<  Specify the number of Arnoldi vectors that are 
                                //!   generated with each iteration (must lie between 
                                //!   2+NEV and N with 2*NEV+1 recommended)
            T* V;               //!<  Memory space to store the Arnoldi vectors (array
                                //!   size N by NCV)
            int LDV;            //!<  Leading dimension of the Arnodli vectors
            int IPARAM[11];     //!<  Sets 11 algorithm parameters:
                                //!  0: Set implicit shifts
                                //!      - 0 user provided
                                //!      - 1 exact shifts from the current Hessenberg matrix
                                //!      - 2 not defined
                                //!  1: Not referenced
                                //!  2: Maximum number of iterations allowed (returns the
                                //!     actual number)
                                //!  3: Block size used in recurrence (only works for 1)
                                //!  4: Specify the number of eigenvalues to include in 
                                //!     the convergence criteria
                                //!  5: Not referenced
                                //!  6: Sets the problem class
                                //!      - 1: Ax = lambda x
                                //!      - 2: Ax = lambda M x (M positive definite)
                                //!      - 3: Ax = lambda M x (M positive semi-definite)
                                //!  7: Holds the number of shifts that the user provides
                                //!     (only used if IPARAM[0] = 2)
                                //!  8: Counts the total number of matrix-vector multiplications
                                //!  9: Counts the total matrix-vector multiplications in the 
                                //!     generalized eigenvalue problem (BMAT = 'G')
                                //! 10: Counts the number of re-orthogonalization steps 
            int IPNTR[14];      //!<  Pointer to mark starting locations for different
                                //!  output types stored in WORKD and WORKL arrays:
                                //!  0: Pointer to current right-vector X in WORKD
                                //!  1: Pointer to current left-vector Y in WORKD
                                //!  2: Pointer to B*X when used in IPARAM 6: 3 mode
                                //!  3: Pointer to next available WORKL location
                                //!  4: Pointer to upper Hessenberg matrix
                                //!  5: Pointer to Ritz array
                                //!  6: Pointer to projected Ritz array
                                //!  7: Pointer to error bounds in WORKL
                                //!  8: Pointer to original Ritz values
                                //!  9: Not used
                                //!  10: Pointer to error bounds
                                //!  11: Pointer to shifts in WORKL
            T* WORKD;           //!<  Working array of dimension 3*N
            T* WORKL;           //!<  Working array of dimension LWORKL
            int LWORKL;         //!<  Dimension of WORKL must be at least  3*NCV**2 + 5*NCV
            double* RWORK;      //!<  Working array of dimension NCV
            int INFO;           //!   Contains output error codes
            //////      A DESCRIPTION OF ADDITIONAL ?neupd_ ARGUMENTS       //////////
            int RECV;           //!<  FALSE: compute Ritz eigenvalues only; TRUE: also
                                //!   compute Ritz vectors
            char HOWMNY;        //!<  'A': compute all NEV Ritz vectors
                                //!   'P': compute all NEV Schur vectors
            int* SELECT;        //!<  Internal work space
            T* D;               //!<  Contains the output Ritz eigenvalues
            T* DR;              //!<  Contains the output Ritz eigenvalues
            T* DI;              //!<  Contains the output Ritz eigenvalues
            T* Z;               //!<  If RECV = true, contains Ritz eigenvectors.
                                //!   Must be of dimension N*NEV
            int LDZ;            //!<  Leading dimension of the Z array 
            T* WORKEV;          //!<  Working space of dimension 2*NCV
            T SIGMA;            //!<  If IPARAM 6: is 3 then contains the "shift value"
            T SIGMAR;           //!<  If IPARAM 6: is 3 then contains the "shift value"
            T SIGMAI;           //!<  If IPARAM 6: is 3 then contains the "shift value"
            //  Remaining arguments are identical to the (p)znaupd_ list    
            //////      ADDITIONAL HIDDEN CHAR LENGTH ARGUMENTS     ////////////////////
            int LEN_BMAT   = 1; //!<  BMAT is a length 1 character
            int LEN_WHICH  = 2; //!<  WHICH is a length 2 character
            int LEN_HOWMNY = 1; //!<  HOWMNY is a length 1 character
            //////      ITERATION COUNTER VARIABLES    //////////
            int iterationCounter = 0;   //!<    Count the number of ARPACK iterations taken
            IDO         = 0;    //  Set value for initial iteration
            BMAT        = m_typeOfProblem;
            N           = nodeDim;
            WHICH[0]    = m_eigenvalueMagnitude;
            WHICH[1]    = m_whatMagnitude;
            NEV         = nbrEigenvalues;
            NCV         = std::min(fullDim, 2*NEV+1+20);
            LDV         = N;
            IPARAM[0]   = 1;        //  Always set ARPACK provided shifts
            IPARAM[2]   = m_maxIterations;
            IPARAM[3]   = m_blockSize;
            IPARAM[4]   = nbrEigenvalues;
            IPARAM[6]   = 1;        //  Always use ARPACK in its standard mode
            if(utilities::is_same<T, dcmplx>::value)
            {
                LWORKL = (3*NCV+5)*NCV;
            }
            else if(utilities::is_same<T, double>::value)
            {
                LWORKL = (3*NCV+6)*NCV;
            }
            INFO        = m_provideInitial;
            RECV        = 1;        //  Compute eigenvectors
            HOWMNY      = 'A';      //  Alternative options not implemented
            LDZ         = N;
            SIGMA       = 0;        //  Set the initial shift value, ARPACK will do the rest
            SIGMAR      = 0;
            SIGMAI      = 0;
            //////      INITIALIZE MEMORY ALLOCATION        ////////////////////////////
            //  Print a warning if we're about to allocate a LOT of working memory
            uint64_t memoryRequirement = (uint64_t)((NEV+1+(1+NCV+3)*N+LWORKL)*sizeof(T));
            if(memoryRequirement>memoryQueryLimit)
            {
                std::cerr<<"\n\tDiagonalizeSymmetricMatrix WARNING: about to allocate "
                <<memoryRequirement/(1024.0*1024.0)<<" MB for znaupd_ routine."
                <<"\n\tPRESS ANY KEY TO CONTINUE."<<std::endl;
                getchar();
            }
            //////      ALLOCATE MEMORY     ////////////////////////////////////////////
            RESID   = new T[N];
            V       = new T[N*NCV];
            WORKD   = new T[3*N];
            WORKL   = new T[LWORKL];
            SELECT  = new int[NCV];
            if(utilities::is_same<T, dcmplx>::value)
            {   
                RWORK = new double[NCV];
                D     = new T[NEV+1];
            }
            else if(utilities::is_same<T, double>::value)
            {
                DR = new T[NEV+1];
                DI = new T[NEV+1];
            }
            Z       = eigenvectors;     //  Set to external memory address
            WORKEV  = new T[2*NCV];
            //////      SET RESID TO THE STARTING VECTOR, IF ONE IS PROVIDED    ////////
            if(m_provideInitial)
            {
                std::ifstream f_initVector;
                f_initVector.open(m_initialVectorFile.c_str(), std::ios::binary);
                if(!f_initVector.is_open())
                {
                    std::cerr<<"\n\tERROR in DiagonalizeSymmetricMatrix: Cannot open initial vector file "<<m_initialVectorFile<<std::endl;
                    std::cerr<<"\n\tUsing random starting vector instead."<<std::endl;
                    INFO = 0;
                }
                else
                {
                    f_initVector.read(reinterpret_cast<char*>(RESID), (long int)N*sizeof(T));
                    f_initVector.close();
                }
            }
            //////      RUN ITERATIVE LANCZOS ALGORITHM       //////////////////////////
            //  Recursively call the (p)znaupd_ routine until the desired eigenvalues
            //  have properly converged
            utilities::cout.AdditionalInfo()<<std::endl;
            do
            {
                utilities::cout.AdditionalInfo()<<"\t ARPACK ITERATION "<<iterationCounter<<"\r";
                fflush(stdout);
                if(utilities::is_same<T, dcmplx>::value)
                {
                    znaupd_(&IDO, &BMAT, &N, &WHICH[0], &NEV, &TOL, RESID, &NCV, V,
                            &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, RWORK, &INFO,
                            &LEN_BMAT, &LEN_WHICH);
                }
                else if(utilities::is_same<T, double>::value)
                {
                    dnaupd_(&IDO, &BMAT, &N, &WHICH[0], &NEV, &TOL, RESID, &NCV, V,
                            &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO,
                            &LEN_BMAT, &LEN_WHICH);
                }
                if(INFO<0)
                {
                    if(utilities::is_same<T, dcmplx>::value)
                    {
                        std::cerr<<"\n\tERROR WITH znaupd_"
                        <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/znaupd.html)"
                        <<std::endl;
                    }
                    else if(utilities::is_same<T, double>::value)
                    {
                        std::cerr<<"\n\tERROR WITH dnaupd_"
                        <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/dnaupd.html)"
                        <<std::endl;
                    }
                    return false;
                }
                else if(1 == INFO)
                {
                    if(utilities::is_same<T, dcmplx>::value)
                    {
                        std::cerr<<"\n\tznaupd_ MESSAGE: RETURNING DUE TO MAXIMUM ITERATIONS REACHED (SET TO "<<m_maxIterations<<")"<<std::endl;
                    }
                    else if(utilities::is_same<T, double>::value)
                    {
                        std::cerr<<"\n\tdnaupd_ MESSAGE: RETURNING DUE TO MAXIMUM ITERATIONS REACHED (SET TO "<<m_maxIterations<<")"<<std::endl;
                    }
                }
                else if(3 == INFO)
                {
                    if(utilities::is_same<T, dcmplx>::value)
                    {
                        std::cerr<<"\n\tznaupd_ MESSAGE: COULD NOT APPLY SHIFTS (TRY INCREASING NCV)"<<std::endl;
                    }
                    else if(utilities::is_same<T, double>::value)
                    {
                        std::cerr<<"\n\tdnaupd_ MESSAGE: COULD NOT APPLY SHIFTS (TRY INCREASING NCV)"<<std::endl;
                    }
                }
                if(1 == IDO || -1 == IDO)
                {
                    //  This output requests a matrix-vector multiplication
                    //  applied to the current Ritz estimates
                    //  Y = matrix * X
                    MatrixVectorMultiply(WORKD+IPNTR[0]-1, WORKD+IPNTR[1]-1);
                }
                ++iterationCounter;
            }
            while(IDO != 99);
            utilities::cout.AdditionalInfo()<<std::endl<<std::endl;
            //  Optionally store the final iteration of the Arnoldi vector
            if(m_storeFinalVector)
            {
                std::ofstream f_finalVector;
                f_finalVector.open(m_finalVectorFile.c_str(), std::ios::binary);
                if(!f_finalVector.is_open())
                {
                    std::cerr<<"\n\tERROR in DiagonalizeSymmetricMatrix: Cannot open final vector file "<<m_finalVectorFile<<std::endl;
                }
                else
                {
                    f_finalVector.write((char*)(V), (long int)N*sizeof(T));
                    f_finalVector.close();
                }
            }
            //////      RUN POST PROCESSING ROUTINE     ///////////////////////////////
            if(utilities::is_same<T, dcmplx>::value)
            {
                zneupd_(&RECV, &HOWMNY, SELECT, D, Z, &LDZ,
                        &SIGMA, WORKEV, &BMAT, &N, &WHICH[0], &NEV, &TOL, RESID,
                        &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL,
                        &LWORKL, RWORK, &INFO, &LEN_HOWMNY, &LEN_BMAT, &LEN_WHICH);
            }
            else if(utilities::is_same<T, double>::value)
            {
                dneupd_(&RECV, &HOWMNY, SELECT, DR, DI, Z, &LDZ,
                        &SIGMAR, &SIGMAI, WORKEV, &BMAT, &N, &WHICH[0], &NEV, &TOL, RESID,
                        &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL,
                        &LWORKL, &INFO, &LEN_HOWMNY, &LEN_BMAT, &LEN_WHICH);
            }
            if(INFO<0)
            {
                if(utilities::is_same<T, dcmplx>::value)
                {
                    std::cerr<<"\n\tERROR WITH zneupd_"
                    <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/zneupd.html)"
                    <<std::endl;
                }
                else if(utilities::is_same<T, double>::value)
                {
                    std::cerr<<"\n\tERROR WITH dneupd_"
                    <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/dneupd.html)"
                    <<std::endl;
                }
                return false;
            }
            for(int i=0; i<nbrEigenvalues; ++i)
            {
                if(utilities::is_same<T, dcmplx>::value)
                {
                    eigenvalues[i] = std::real(D[i]);
                }
                else if(utilities::is_same<T, double>::value)
                {
                    eigenvalues[i] = std::real(DR[i]);
                }
            }
            //////      CLEAR UP MEMORY     ////////////////////////////////////////////
            delete[] RESID;
            delete[] V;
            delete[] WORKD;
            delete[] WORKL;
            delete[] SELECT;
            if(utilities::is_same<T, dcmplx>::value)
            {
                delete[] RWORK;
                delete[] D;
            }
            else if(utilities::is_same<T, double>::value)
            {
                delete[] DR;
                delete[] DI;
            }
            delete[] WORKEV;
            //////      PRINT SUMMARY OF CALCULATION DATA    ///////////////////////////
            utilities::cout.SecondaryOutput()<<"\n\tARPACK NCV:\t\t\t"<<NCV<<std::endl;
            utilities::cout.SecondaryOutput()<<"\tARPACK ITERATIONS:\t\t"<<iterationCounter<<std::endl;
            utilities::cout.SecondaryOutput()<<"\tARPACK NO. MATRIX-VECTOR CALLS:\t"<<iterationCounter<<std::endl;
            return true;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief This function calls an PARPACK diagonalization routine to obtain a
        //!	approximate eigenvalues/eigenvectors of a complex hermitian matrix
        //! or real symmetric matrix
        //!
        //! Returns selected Ritz values/vectors (approximate eigenvalues/eigenvectors) 
        //! to the memory address specified by the eigenvalues/eigenvectors variables. 
        //!
        //! NOTE: The eigenvalues will be reproduced on each node, but the eigenvectors
        //! will be distributed over nodes (the first x rows are stored on node 1
        //! , then next x rows on node 2 etc. where x = dimension /# nodes)
        //!
        //! NOTE: shift invert mode is not implemented in parallel, since I have
        //! not incorperated a parallel Cholesky factorization and linear solver
        //! method
        ////////////////////////////////////////////////////////////////////////////////
        template <typename T>
        bool DiagonalizeSymmetricMatrix(
        const std::function<void(T* inVector, T* outVector, MpiWrapper& mpi)>& MatrixVectorMultiply,
                                    //!<    Function to perform matrix-vector operation
        T* eigenvectors,            //!<    A memory address (N by M in size) to store
                                    //!     the resulting eigenvectors
        double* eigenvalues,        //!<    The memory address (N in size) to store
                                    //!     the resulting eigenvalues
        const int nodeDim,          //!<    Dimension of the matrix/eigenvectors stored on each node
        const int fullDim,          //!<    Full dimension of the matrix
        const int nbrEigenvalues,   //!<    Number of eigenvalues requested
        MpiWrapper& mpi)            //!<    Instance of the MPI wrapper class
        {
            //////      A DESCRIPTION OF ALL p?naupd_ ARGUMENTS     ///////////////////
            int COMM;           //!<  MPI communicator
            int IDO;            //!<  Must be set to zero initially,
                                //!   afterwards specifies the type of operation performed
                                //!   and value 99 tells us the iterative procedure has finished
            char BMAT;          //!<  'I' sets a standard eigenvalue problem, 'G' sets a 
                                //!  generalized eigenvalue problem instead
            int N;              //!<  Dimension of the eigenvalue problem
            char WHICH[2];      //!<  Option for range of eigenvalue returned:
                                //!  'LM' -> eigenvalues of largest magnitude
                                //!  'SM' -> eigenvalues of smallest magnitude
                                //!  'LR' -> eigenvalues of largest real part
                                //!  'SR' -> eigenvalues of smallest real part
                                //!  'LI' -> eigenvalues of largest imaginary part
                                //!  'SI' -> eigenvalues of smallest imaginary part
            int NEV;            //!<   Number of eigenvalues/eigenvectors to find
            double TOL;         //!<  Stopping tolerance (defaults to machine precision)
            T* RESID;           //!<  A dimension-N vector containing the residuals
            int NCV;            //!<  Specify the number of Arnoldi vectors that are 
                                //!   generated with each iteration (must lie between 
                                //!   2+NEV and N with 2*NEV+1 recommended)
            T* V;               //!<  Memory space to store the Arnoldi vectors (array
                                //!   size N by NCV)
            int LDV;            //!<  Leading dimension of the Arnodli vectors
            int IPARAM[11];     //!<  Sets 11 algorithm parameters:
                                //!  0: Set implicit shifts
                                //!      - 0 user provided
                                //!      - 1 exact shifts from the current Hessenberg matrix
                                //!      - 2 not defined
                                //!  1: Not referenced
                                //!  2: Maximum number of iterations allowed (returns the
                                //!     actual number)
                                //!  3: Block size used in recurrence (only works for 1)
                                //!  4: Specify the number of eigenvalues to include in 
                                //!     the convergence criteria
                                //!  5: Not referenced
                                //!  6: Sets the problem class
                                //!      - 1: Ax = lambda x
                                //!      - 2: Ax = lambda M x (M positive definite)
                                //!      - 3: Ax = lambda M x (M positive semi-definite)
                                //!  7: Holds the number of shifts that the user provides
                                //!     (only used if IPARAM[0] = 2)
                                //!  8: Counts the total number of matrix-vector multiplications
                                //!  9: Counts the total matrix-vector multiplications in the 
                                //!     generalized eigenvalue problem (BMAT = 'G')
                                //! 10: Counts the number of re-orthogonalization steps 
            int IPNTR[14];      //!<  Pointer to mark starting locations for different
                                //!  output types stored in WORKD and WORKL arrays:
                                //!  0: Pointer to current right-vector X in WORKD
                                //!  1: Pointer to current left-vector Y in WORKD
                                //!  2: Pointer to B*X when used in IPARAM 6: 3 mode
                                //!  3: Pointer to next available WORKL location
                                //!  4: Pointer to upper Hessenberg matrix
                                //!  5: Pointer to Ritz array
                                //!  6: Pointer to projected Ritz array
                                //!  7: Pointer to error bounds in WORKL
                                //!  8: Pointer to original Ritz values
                                //!  9: Not used
                                //!  10: Pointer to error bounds
                                //!  11: Pointer to shifts in WORKL
            T* WORKD;           //!<  Working array of dimension 3*N
            T* WORKL;           //!<  Working array of dimension LWORKL
            int LWORKL;         //!<  Dimension of WORKL must be at least  3*NCV**2 + 5*NCV
            double* RWORK;      //!<  Working array of dimension NCV
            int INFO;           //!   Contains output error codes
            //////      A DESCRIPTION OF ADDITIONAL p?neupd_ ARGUMENTS       //////////
            int RECV;           //!<  FALSE: compute Ritz eigenvalues only; TRUE: also
                                //!   compute Ritz vectors
            char HOWMNY;        //!<  'A': compute all NEV Ritz vectors
                                //!   'P': compute all NEV Schur vectors
            int* SELECT;        //!<  Internal work space
            T* D;               //!<  Contains the output Ritz eigenvalues
            T* DR;              //!<  Contains the output Ritz eigenvalues
            T* DI;              //!<  Contains the output Ritz eigenvalues
            T* Z;               //!<  If RECV = true, contains Ritz eigenvectors.
                                //!   Must be of dimension N*NEV
            int LDZ;            //!<  Leading dimension of the Z array
            T* WORKEV;          //!<  Working space of dimension 2*NCV
            T SIGMA;            //!<  If IPARAM 6: is 3 then contains the "shift value"
            T SIGMAR;           //!<  If IPARAM 6: is 3 then contains the "shift value"
            T SIGMAI;           //!<  If IPARAM 6: is 3 then contains the "shift value"
            //////      ADDITIONAL HIDDEN CHAR LENGTH ARGUMENTS     ////////////////////
            int LEN_BMAT   = 1; //!<  BMAT is a length 1 character
            int LEN_WHICH  = 2; //!<  WHICH is a length 2 character
            int LEN_HOWMNY = 1; //!<  HOWMNY is a length 1 character
            //////      ITERATION COUNTER VARIABLES    //////////
            int iterationCounter = 0;   //!<    Count the number of ARPACK iterations taken
            COMM        = MPI_Comm_c2f(mpi.m_comm); //!<  Converted MPI comm handle
            IDO         = 0;    //!<  Set value for initial iteration
            BMAT        = m_typeOfProblem;
            N           = nodeDim;
            WHICH[0]    = m_eigenvalueMagnitude;
            WHICH[1]    = m_whatMagnitude;
            NEV         = nbrEigenvalues;
            NCV         = std::min(fullDim, 2*NEV+1+20);
            LDV         = N;
            IPARAM[0]   = 1;    //!<  Always set ARPACK provided shifts
            IPARAM[2]   = m_maxIterations;
            IPARAM[3]   = m_blockSize;
            IPARAM[4]   = nbrEigenvalues;
            IPARAM[6]   = 1;    //!<  Always use ARPACK in its standard mode
            if(utilities::is_same<T, dcmplx>::value)
            {
                LWORKL = (3*NCV+5)*NCV;
            }
            else if(utilities::is_same<T, double>::value)
            {
                LWORKL = (3*NCV+6)*NCV;
            }
            INFO        = m_provideInitial;
            RECV        = 1;        //  Compute eigenvectors
            HOWMNY      = 'A';      //  Alternative options not implemented
            LDZ         = N;
            SIGMA       = 0;        //  Set the initial shift value, ARPACK will do the rest
            SIGMAR      = 0;
            SIGMAI      = 0;
            //////      INITIALIZE MEMORY ALLOCATION        ////////////////////////////
            uint64_t memoryRequirement = (uint64_t)((NEV+1+(1+NCV+3)*N+LWORKL)*sizeof(T));
            if(memoryRequirement>memoryQueryLimit)
            {
                std::cerr<<"\n\tDiagonalizeSymmetricMatrix WARNING: about to allocate "
                <<memoryRequirement/(1024.0*1024.0)<<" MB for znaupd_ routine."
                <<"\n\tPRESS ANY KEY TO CONTINUE."<<std::endl;
                getchar();
            }
            MPI_Barrier(mpi.m_comm);
            //////      ALLOCATE MEMORY     ////////////////////////////////////////////
            RESID   = new T[N];
            V       = new T[N*NCV];
            WORKD   = new T[3*N];
            WORKL   = new T[LWORKL];
            SELECT  = new int[NCV];
            if(utilities::is_same<T, dcmplx>::value)
            {   
                RWORK = new double[NCV];
                D     = new T[NEV+1];
            }
            else if(utilities::is_same<T, double>::value)
            {
                DR = new T[NEV+1];
                DI = new T[NEV+1];
            }
            Z       = eigenvectors;     //  Set to external memory address
            WORKEV  = new T[2*NCV];
            //////      SET RESID TO THE STARTING VECTOR, IF ONE IS PROVIDED    ////////
            if(m_provideInitial)
            {
                //  Open file containing initial vector 
                //  (when run in parallel, the file stores the part of the 
                //  initial vector on each node)
                std::ifstream f_initVector;
                f_initVector.open(m_initialVectorFile.c_str(), std::ios::binary);
                if(!f_initVector.is_open())
                {
                    if(0 == mpi.m_id)
                    {
                        std::cerr<<"\n\tERROR in DiagonalizeSymmetricMatrix: Cannot open initial vector file "<<m_initialVectorFile<<std::endl;
                        std::cerr<<"\n\tUsing random starting vector instead."<<std::endl;
                    }
                    INFO = 0;
                }
                else
                {
                    f_initVector.read(reinterpret_cast<char*>(RESID),(long int)N*sizeof(T));
                    f_initVector.close();
                }
            }
            //////      RUN ITERATIVE LANCZOS ALGORITHM       //////////////////////////
            //  Recursively call the (p)znaupd_ routine until the desired eigenvalues
            //  have properly converged
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<std::endl;
            }
            do
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\t ARPACK ITERATION "<<iterationCounter<<"\r";
                    fflush(stdout);
                }
                if(utilities::is_same<T, dcmplx>::value)
                {
                    pznaupd_(&COMM, &IDO, &BMAT, &N, &WHICH[0], &NEV, &TOL, RESID,
                             &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, RWORK, &INFO,
                             &LEN_BMAT, &LEN_WHICH);
                }
                else if(utilities::is_same<T, double>::value)
                {
                    pdnaupd_(&COMM, &IDO, &BMAT, &N, &WHICH[0], &NEV, &TOL, RESID,
                             &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO,
                             &LEN_BMAT, &LEN_WHICH);
                }
                if(INFO<0)
                {
                    if(utilities::is_same<T,dcmplx>::value)
                    {
                        std::cerr<<"\n\tERROR WITH pznaupd_"
                        <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/pznaupd.html)"
                        <<std::endl;
                    }
                    else if(utilities::is_same<T,double>::value)
                    {
                        std::cerr<<"\n\tERROR WITH pdnaupd_"
                        <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/pdnaupd.html)"
                        <<std::endl;
                    }

                    return false;
                }
                else if(1 == INFO)
                {
                    if(utilities::is_same<T,dcmplx>::value)
                    {
                        std::cerr<<"\n\tpznaupd_ MESSAGE: RETURNING DUE TO MAXIMUM ITERATIONS REACHED (SET TO "<<m_maxIterations<<")"<<std::endl;
                    }
                    else if(utilities::is_same<T,double>::value)
                    {
                        std::cerr<<"\n\tpdnaupd_ MESSAGE: RETURNING DUE TO MAXIMUM ITERATIONS REACHED (SET TO "<<m_maxIterations<<")"<<std::endl;
                    }
                }
                else if(3 == INFO)
                {
                    if(utilities::is_same<T,dcmplx>::value)
                    {
                        std::cerr<<"\n\tpznaupd_ MESSAGE: COULD NOT APPLY SHIFTS (TRY INCREASING NCV)"<<std::endl;
                    }
                    else if(utilities::is_same<T,double>::value)
                    {
                        std::cerr<<"\n\tpdnaupd_ MESSAGE: COULD NOT APPLY SHIFTS (TRY INCREASING NCV)"<<std::endl;
                    }
                }
                if(1 == IDO || -1 == IDO)
                {
                    //  This output requests a matrix-vector multiplication
                    //  applied to the current Ritz estimates
                    //  Y = matrix * X
                    MatrixVectorMultiply(WORKD+IPNTR[0]-1, WORKD+IPNTR[1]-1, mpi);
                }
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    ++iterationCounter;
                }
            }
            while(IDO != 99);
            MPI_Barrier(mpi.m_comm);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<std::endl<<std::endl;
            }
            //  Optionally store the final iteration of the Arnoldi vector
            if(m_storeFinalVector)
            {
                std::ofstream f_finalVector;
                f_finalVector.open(m_finalVectorFile.c_str(), std::ios::binary);
                if(!f_finalVector.is_open())
                {
                    std::cerr<<"\n\tERROR in DiagonalizeSymmetricMatrix: Cannot open final vector file "<<m_finalVectorFile<<std::endl;
                }
                else
                {
                    f_finalVector.write((char*)(V), (long int)N*sizeof(T));
                    f_finalVector.close();
                }
            }
            //////      RUN POST PROCESSING ROUTINE     ///////////////////////////////
            if(utilities::is_same<T, dcmplx>::value)
            {
                pzneupd_(&COMM, &RECV, &HOWMNY, SELECT, D, Z, &LDZ,
                         &SIGMA, WORKEV, &BMAT, &N, &WHICH[0], &NEV,
                         &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR,
                         WORKD, WORKL, &LWORKL, RWORK, &INFO, &LEN_HOWMNY, &LEN_BMAT, &LEN_WHICH);
            }
            else if(utilities::is_same<T, double>::value)
            {
                pdneupd_(&COMM, &RECV, &HOWMNY, SELECT, DR, DI, Z, &LDZ,
                         &SIGMAR, &SIGMAI, WORKEV, &BMAT, &N, &WHICH[0], &NEV,
                         &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR,
                         WORKD, WORKL, &LWORKL, &INFO, &LEN_HOWMNY, &LEN_BMAT, &LEN_WHICH);
            }   
            if(INFO<0)
            {
                if(utilities::is_same<T, dcmplx>::value)
                {
                    std::cerr<<"\n\tERROR WITH pzneupd_"
                    <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/pzneupd.html)"
                    <<std::endl;
                }
                else if(utilities::is_same<T, double>::value)
                {
                    std::cerr<<"\n\tERROR WITH pdneupd_"
                    <<"\n\t ERROR CODE WAS "<<INFO<<" (see http://www.mathkeisan.com/usersguide/man/pdneupd.html)"
                    <<std::endl;
                }
                return false;
            }
            for(int i=0; i<nbrEigenvalues; ++i)
            {    
                if(utilities::is_same<T, dcmplx>::value)
                {
                    eigenvalues[i] = std::real(D[i]);
                }
                else if(utilities::is_same<T, double>::value)
                {
                    eigenvalues[i] = std::real(DR[i]);
                }
            }
            //////      CLEAR UP MEMORY     ////////////////////////////////////////////
            delete[] RESID;
            delete[] V;
            delete[] WORKD;
            delete[] WORKL;
            delete[] SELECT;
            if(utilities::is_same<T, dcmplx>::value)
            {
                delete[] RWORK;
                delete[] D;
            }
            else if(utilities::is_same<T, double>::value)
            {
                delete[] DR;
                delete[] DI;
            }
            delete[] WORKEV;
            //////      PRINT SUMMARY OF CALCULATION DATA    ///////////////////////////
            if(0 == mpi.m_id)    //  FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\tPARPACK NCV:\t\t\t"<<NCV<<std::endl;
                utilities::cout.SecondaryOutput()<<"\tPARPACK ITERATIONS:\t\t"<<iterationCounter<<std::endl;
                utilities::cout.SecondaryOutput()<<"\tPARPACK NO. MATRIX-VECTOR CALLS:\t"<<iterationCounter<<std::endl;
            }
            return true;
        }
    };  //  End class ArpackWrapper
    }   //  End namespace linearAlgebra 
}   //  End namespace utilities
#endif
