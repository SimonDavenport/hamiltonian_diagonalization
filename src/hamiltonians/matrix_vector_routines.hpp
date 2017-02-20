////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 26/06/2015
//!                                                                             
//!	 \file
//!     This file contains a class to implement a variety of different 
//!     matrix-vector routines, to be used in Lancoz type algorithms
//!
//!                    Copyright (C) 2015 Simon C Davenport
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

#ifndef _MATRIX_VECTOR_ROUTINES_HPP_INCLUDED_
#define _MATRIX_VECTOR_ROUTINES_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

//  Functions to calculate the Hamiltonian
#include "../hamiltonians/spinless_fermion_hamiltonian.hpp"
#include "../hamiltonians/two_level_spinless_fermion_hamiltonian.hpp"

//  Include various utility classes and definitions
#include "../utilities/wrappers/mpi_wrapper.hpp"

#include <functional>                       //  For passing matrix-vector routine

#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{

template <typename T,class L1,class L2>
class MatrixVectorFunction
{
    friend class HamiltonianBase<T,FermionFockBasis>;
    friend class HamiltonianBase<T,TwoLevelFermionFockBasis>;

    private:

    std::vector<int> dataDistribution;      //!<    Array containing a list of 
                                            //!     vector dimensions on each node
    std::vector<MPI_Comm> commGroups;       //!<    A list of minimal MPI communicator
                                            //!     groups
    std::vector<T> inputVectorBuffer;       //!<    Buffer used in input vector 
                                            //!     MPI communication
    std::vector<T> outputVectorBuffer;      //!<    Buffer used in output vector
                                            //!     MPI communication
    char UPLO;                              //!<    Flag to state what part of the 
                                            //!     Hermitian matrix is stored
    int matrixPower;                        //!<    Matrix power in preconditioner
    std::vector<double> a,b;                //!<    Shift and stretch parameters
                                            //!     in preconditioner  
    utilities::CrsSparseMatrix<T>* matrix;  //!<    Optional pointer to the CRS
                                            //!     sparse matrix
    std::vector<T> tempVector;              //!<    Temporary storage to be used
                                            //!     if matrixPower>1

    SpinlessFermionHamiltonian<T>* hamiltonian;
                                            //!<    Optional pointer to a spinless
                                            //!     Hamiltonian class
    
    TwoLevelSpinlessFermionHamiltonian<T>* twoLevelHamiltonian;
                                            //!<    Optional pointer to a two-level
                                            //!     spinless Hamiltonian class
    public:

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor from a spinless fermion Hamiltonian
    //! 
    ////////////////////////////////////////////////////////////////////////////////
    
    MatrixVectorFunction(
        SpinlessFermionHamiltonian<T>& hamiltonian,
                                            //!<    Instance of the Hamiltonian class
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        this->UPLO = hamiltonian.m_arpackWrapper.GetUpperOrLower();

        hamiltonian.m_arpackWrapper.GetPreconditioner(this->matrixPower,this->a,this->b);
        
        if(this->matrixPower>1)
        {
            tempVector.resize(hamiltonian.m_data.m_nodeDim);
        }
        
        this->matrix = &(hamiltonian.m_crsSparseMatrix);
        
        this->hamiltonian = &hamiltonian;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor from a two level spinless fermion Hamiltonian
    //! 
    ////////////////////////////////////////////////////////////////////////////////
    
    MatrixVectorFunction(
        TwoLevelSpinlessFermionHamiltonian<T>& hamiltonian,
                                            //!<    Instance of the Hamiltonian class
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        this->UPLO = hamiltonian.m_arpackWrapper.GetUpperOrLower();

        hamiltonian.m_arpackWrapper.GetPreconditioner(this->matrixPower,this->a,this->b);
        
        //if(this->matrixPower>1)
        {
            tempVector.resize(hamiltonian.m_data.m_nodeDim);
        }
        
        this->matrix = &(hamiltonian.m_crsSparseMatrix);
        
        this->twoLevelHamiltonian = &hamiltonian;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Perform a simple matrix-vector multiplication routine.
    //! 
    ////////////////////////////////////////////////////////////////////////////////
    
    void Simple(
        T* x,                               //!<    Input vector on the current node
        T* y,                               //!<    Output vector on the current node
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        utilities::linearAlgebra::ParallelSymmetricMatrixVectorMultiply(x,y,&(this->dataDistribution),this->matrix,&(this->UPLO),(T)this->a[0],(T)this->b[0],&(this->commGroups),this->inputVectorBuffer.data(),this->outputVectorBuffer.data(),this->inputVectorBuffer.size(),mpi);
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Simulate obtaining the eigenvalues of a matrix to some power
    //! by repeatedly appling the matrix-vector routine
    //! 
    ////////////////////////////////////////////////////////////////////////////////
    
    void Preconditioned(
        T* x,                               //!<    Input vector on the current node
        T* y,                               //!<    Output vector on the current node
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        this->Simple(x,y,mpi);
    
        for(int i=1;i<matrixPower;i++)
        {
            memcpy(tempVector.data(),y,tempVector.size()*sizeof(T));

            utilities::linearAlgebra::ParallelSymmetricMatrixVectorMultiply(tempVector.data(),y,&(this->dataDistribution),this->matrix,&(this->UPLO),(T)this->a[0],(T)this->b[0],&(this->commGroups),this->inputVectorBuffer.data(),this->outputVectorBuffer.data(),this->inputVectorBuffer.size(),mpi);
        }
    }
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief matrix-vector function involving recalculation of the matrix at
    //! each step
    //!
    //! TODO debug this function - it does not currently work!
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void Recalculate(
        L1* quadraticLookUpTables,          //!<    Class container for look-up tables
        L2* quarticLookUpTables,            //!<    Class container for look-up tables
        T* x,                               //!<    Input vector on the current node
        T* y,                               //!<    Output vector on the current node
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        fock_t nbrRowsPerSweep = 2;
        fock_t nbrSweeps = std::ceil((double)this->hamiltonian->m_data.m_nodeDim/nbrRowsPerSweep);
        fock_t currRow = 0;
        fock_t nbrNonZeros = 10000;
        
        //  Generate and store in CRS format a subset of all the matrix rows then
        //  perform a reduced matrix-vector multiplication for that subset of
        //  rows and store the outcome in a buffer

        this->hamiltonian->m_crsSparseMatrix.Clear();
        
        this->hamiltonian->m_crsSparseMatrix.Initialize(this->hamiltonian->m_data.m_nodeDim,this->hamiltonian->m_data.m_fockSpaceDim,nbrNonZeros);

        for(fock_t sweep=0;sweep<nbrSweeps;++sweep,currRow+=nbrRowsPerSweep)
        {
            this->hamiltonian->Add_CdC_Terms(quadraticLookUpTables,currRow,currRow+nbrRowsPerSweep,mpi);
            this->hamiltonian->Add_CdCdCC_Terms(quarticLookUpTables,currRow,currRow+nbrRowsPerSweep,mpi);
            
            this->Simple(x,y,mpi);
            
            for(fock_t i=0;i<this->hamiltonian->m_data.m_nodeDim;++i)
            {
                tempVector[i] += y[i];
            }
        }
        
        //  Copy the accumulated vector into the output buffer
        
        memcpy(y,tempVector.data(),tempVector.size()*sizeof(T));
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
};

}   //  End diagonalization namespace

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#endif
