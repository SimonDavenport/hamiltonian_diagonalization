////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a class to generate a bitwise encoded occupation
//!     basis for spinless Fermion orbitals and generate a sparse Hermitian
//!     Hamiltonian optionally containing both quadratic and quartic interactions. 
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

#ifndef _SPINLESS_FERMION_HAMILTONIAN_HPP_INCLUDED_
#define _SPINLESS_FERMION_HAMILTONIAN_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "hamiltonian_base.hpp"
#include "../utilities/mathematics/binary_number_tools.hpp"
#include "../utilities/mathematics/fermions.hpp"
#include "../basis/fermion_fock_basis.hpp"    
#if _DEBUG_
#include "../utilities/general/debug.hpp"
#include <bitset>
#endif

namespace diagonalization
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief The SpinlessFermionHamiltonian class generates and contains a 
    //! set of basis states that are a bitwise representation of fermions
    //! occupying a set of orbitals
    //!
    //!	Fock states are represented in ascending order c^+_0 c^+_2 c^+_5 etc. 
    //! and are encoded in uint64_t binaries
    //!
    //! The template parameter T denotes the data type of the matrix elements
    ////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    class SpinlessFermionHamiltonian : public HamiltonianBase<T, FermionFockBasis>
    {
	    public:
        //!
        //! Initialization function to override default construction
        //! 
        void Initialize(
            const iSize_t nbrParticles,             //!<    Number of particles
            const iSize_t nbrOrbitals,              //!<    Number of orbitals
            const utilities::storageMethod_t storageMethod, 
                                                    //!<    Matrix storage method
            utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
        {
            this->m_data = HamiltonianData(nbrParticles, nbrOrbitals, storageMethod);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                if(0 >= this->m_data.m_nbrParticles)
                {
                    std::cerr<<"\n\tERROR: MUST SPECIFY nbr particles > 0"<<std::endl;
                    mpi.m_exitFlag=true;
                }
                if(0 >= this->m_data.m_nbrOrbitals)
                {
                    std::cerr<<"\n\tERROR: MUST SPECIFY nbr orbitals > 0"<<std::endl;
                    mpi.m_exitFlag=true;
                }
                if(this->m_data.m_nbrParticles >= 64)
                {
                    std::cerr<<"\n\n\tERROR: CANNOT ADDRESS 64 PARTICLES OR MORE.\n"<<std::endl;
                    mpi.m_exitFlag=true;
                }
                if(this->m_data.m_nbrParticles > this->m_data.m_nbrOrbitals)
                {
	                std::cerr<<"\n\n\tERROR: NBR PARTICLES > NBR ORBITALS - NOT ALLOWED\n"<<std::endl;
	                mpi.m_exitFlag=true;
                }
            }
            mpi.ExitFlagTest();
            this->m_data.SetDataDimensions(this->m_fockBasis.GetReducedFockSpaceDimension(
                                           this->m_data.m_nbrParticles, this->m_data.m_nbrOrbitals, mpi),
                                           mpi.m_nbrProcs, mpi);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            { 
                utilities::cout.SecondaryOutput()<<"\n\tHAMILTONIAN INITIALIZED WITH THE FOLLOWING PARAMETERS:\n"<<std::endl;
                utilities::cout.SecondaryOutput()<<"\t\tNBR PARTICLES:\t\t"<<this->m_data.m_nbrParticles<<std::endl;
                utilities::cout.SecondaryOutput()<<"\t\tNBR ORBITALS:\t\t"<<this->m_data.m_nbrOrbitals<<std::endl;
                this->BaseInitialize(mpi);
	        }
	        //  Set estimated number of non zeros to be the twice fock space dimension
	        //  on a given node
	        this->AllocateMatrixMemory(this->m_data.m_storageMethod, 2*this->m_data.m_nodeDim, _PARALLEL_, mpi);
	        mpi.ExitFlagTest();
        }
        
        //!
        //! Clear all class data structures
        //!
        void Clear()
        {
            this->BaseClear();
        }

        //!
        //! Destructor
        //!
        ~SpinlessFermionHamiltonian()
        {}
 
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A template function to add all valid quartic terms of the form
        //!  E_{k1,k2} c_k1^+ c_k2
        //!
        //! The function argument must be a class containing the following three methods:
        //!
        //! - GetMaxKCount()
        //!
        //! Return the maximum number of k1 values that will be returned given any k2
        //! due to any kind of conservation law
        //!
        //! - GetK1(kRetrieveBuffer,nbrK1,k2)
        //!
        //! To return a list of k1 values related to a given k2 by some conservation
        //! law (where each k represents the pair kx,ky on a 2D grid).
        //!
        //! - GetEkk(k1,k2)
        //!
        //! To return a single Ekk coefficient for the specified k indices
        //!
        //! This function MUST be called on all nodes if run in parallel
        ////////////////////////////////////////////////////////////////////////////////
        template <class L>
        void Add_CdC_Terms(
            L* termTables,              //!<    Class container for terms tables
            fock_t minRow,              //!<    Minimum row index to calculate for
            fock_t maxRow,              //!<    Maximum row index to calculate for
            utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
        {
            if(this->m_data.m_matrixAllocated)
	        {
	            utilities::Timer timer;
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\n\t- ADDING C^+ C TERMS...\n"<<std::endl;
                    timer.Start();     //  Time the operation
                }
                //  For parallel processes, divide up tasks amongst available nodes
                mpi.DivideTasks(mpi.m_id, this->m_data.m_fockSpaceDim,
                                mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
                fock_t lowerRowLimit = std::max((fock_t)mpi.m_firstTask, mpi.m_firstTask+minRow);
                fock_t upperRowLimit = std::min((fock_t)mpi.m_lastTask, mpi.m_firstTask+maxRow);
                //  Reserve a memory block with which to retrieve the momentum conserving list of k-values
                kState_t kRetrieveBuffer[termTables->GetMaxKCount()];
                //  Determine the starting Fock state on each node
                std::vector<fock_t>::const_iterator it_state;
                fock_t currState = this->m_fockBasis.GetFirstState(lowerRowLimit, this->m_data.m_nbrParticles, it_state);
                utilities::LoadBar progress;
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    progress.Initialize(this->m_data.m_nodeDim);
                }
                for(int i=lowerRowLimit; i<=upperRowLimit; ++i)
                {
                    fock_t inState = currState;
                    //  First, we can restrict the sum to k2 values that do no annihilate
                    //  the inState
                    fock_t tempK2 = inState;
                    while(tempK2!=0)    //    N iterations for N particles
                    {
                        fock_t k2Occupied = tempK2 & -tempK2; 
                        kState_t k2 = utilities::binary::HammingWeight64(k2Occupied-1);
                        iSize_t nbrK1;
                        termTables->GetK1(kRetrieveBuffer, nbrK1, k2);             
                        //  Loop over the possible k1 values
                        for(iSize_t j=0; j<nbrK1; ++j)
                        {
                            kState_t k1 = kRetrieveBuffer[j];
                            T matrixElement = 0.0;
                            fock_t outIndex = i;
                            if(k1==k2)    //  Case of diagonal elements
                            {
                                matrixElement = termTables->GetEkk(k1, k2);
                            }
                            else
                            {
                                //  Zero if k1 occurs in the inState
                                fock_t k1Occupied = (fock_t)1 << k1;
                                if(k1Occupied & ~inState)
                                {
                                    //  Act on the current in state with CdC
                                    int sign;
                                    fock_t outState = utilities::Evaluate_CdC(sign, k1Occupied, k2Occupied, inState,
                                                                              this->m_data.m_highestOrbital);
                                    if(0 != sign && outState>=inState)
                                    {
                                        //  If non zero then calculate the matrix element
                                        matrixElement = termTables->GetEkk(k1, k2)*(double)sign;
                                        outIndex = this->m_fockBasis.FindFockStateIndex(mpi.m_firstTask, outState);
                                    }
                                }
                            }
                            //  Don't bother to add the matrix element if it is zero
                            if(outIndex < this->m_data.m_fockSpaceDim)
                            {
                                if((utilities::is_same<dcmplx,T>::value && abs(matrixElement)>_SPARSE_ZERO_TOL_) || (utilities::is_same<double, T>::value && fabs(matrixElement)>_SPARSE_ZERO_TOL_))
                                {
                                    if(utilities::_DENSE_ == this->m_data.m_storageMethod)
                                    {
                                        this->m_fullMatrix[(i-mpi.m_firstTask)*this->m_data.m_fockSpaceDim+outIndex] += matrixElement;
                                    }  
                                    else if(utilities::_SPARSE_MAPPED_ == this->m_data.m_storageMethod)
                                    {
                                        this->m_sparseMap.Insert((i-mpi.m_firstTask),outIndex) += matrixElement;
                                    }
                                    else if(utilities::_SPARSE_CRS_ == this->m_data.m_storageMethod)
                                    {
                                        this->m_crsSparseMatrix.Insert((i-mpi.m_firstTask),outIndex)  += matrixElement;
                                    }
                                }
                            }
                        }
                        tempK2 = tempK2 ^ k2Occupied;     //  Remove right-most bit and proceed to the next iteration
                    }
                    currState = this->m_fockBasis.GetNextState(currState, it_state);
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {
                        progress.Display(i+1);
                    }
                }
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<std::endl;
                    timer.Stop();
                }
            }
            else
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::cerr<<"\n\tERROR in Add_CdC_Terms: Matrix memory not allocated!"<<std::endl;
                    mpi.m_exitFlag = true;
                }
            }
            mpi.ExitFlagTest();
	        return;
        }
  
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A template function to add all valid quartic terms of the form 
        //! V_{k1,k2,k3,k4} c_k1^+ c_k2^+ c_k3 c_k4
        //!
        //! The function argument must be a class containing the following three methods:
        //!
        //! - GetMaxKCount()
        //!
        //! Return the maximum number of k1 values that will be returned given any 
        //! k2,k3,k4 due to any kind of conservation law
        //!
        //! - GetK1(kRetrieveBuffer,nbrK1,k2,k3,k4)
        //!
        //! To return a list of k1 values related to a given k2,k3,k4 by some 
        //! conservation law (where each k represents the pair kx,ky on a 2D grid).
        //!
        //! - GetVkkkk(k1,k2,k3,k4)
        //!
        //! To return a single Vkkkk coefficient for the specified k indices
        //!
        //! This function MUST be called on all nodes if run in parallel
        ////////////////////////////////////////////////////////////////////////////////
        template <class L>
        void Add_CdCdCC_Terms(
            L* termTables,              //!<    Class container for terms tables
            fock_t minRow,              //!<    Minimum row index to calculate for
            fock_t maxRow,              //!<    Maximum row index to calculate for
            utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
        {
	        if(this->m_data.m_matrixAllocated)
	        {
	            utilities::Timer timer;
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\n\t- ADDING C^+ C^+ C C TERMS..."<<std::endl;

                    timer.Start();     //  Time the operation
                }
                mpi.DivideTasks(mpi.m_id, this->m_data.m_fockSpaceDim, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
                fock_t lowerRowLimit = std::max((fock_t)mpi.m_firstTask, mpi.m_firstTask+minRow);
                fock_t upperRowLimit = std::min((fock_t)mpi.m_lastTask, mpi.m_firstTask+maxRow);
                //  Reserve a memory block with which to retrieve the momentum conserving list of k-values
                kState_t kRetrieveBuffer[termTables->GetMaxKCount()];
                std::vector<fock_t>::const_iterator it_state;
                fock_t currState = this->m_fockBasis.GetFirstState(lowerRowLimit, this->m_data.m_nbrParticles, it_state);
                utilities::LoadBar progress;
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    progress.Initialize(this->m_data.m_nodeDim);
                }
                //  Calculate the sum over Vkkkk terms
                for(int i=lowerRowLimit; i<=upperRowLimit; ++i)
                {
                    fock_t inState = currState;
                    //  Determine the k1,k2,k3,k4 orbitals that lead to a non-zero contribution, 
                    //  First, we can restrict the sum to k3,k4 values that do no annihilate
                    //  the inState and that satisfy k3!=k4
                    fock_t tempK4 = inState;
                    const fock_t mask = ((fock_t)1 << this->m_data.m_nbrOrbitals)-1;
                    while(tempK4!=0)    //    N iterations for N particles
                    {
                        fock_t k4Occupied = tempK4 & -tempK4;      //  Isolate the right-most bit
                        kState_t k4 = utilities::binary::HammingWeight64(k4Occupied-1);
                        //  Now process the in state with the k4 orbital removed
                        fock_t tempK3 = inState ^ k4Occupied;
                        while(tempK3!=0)    //    N iterations for N particles
                        {
                            fock_t k3Occupied = tempK3 & -tempK3;      //  Isolate the right-most bit
                            kState_t k3 = utilities::binary::HammingWeight64(k3Occupied-1);
                            if(k3 != k4)        //  For fermions, a squared creation operator gives zero
                            {
                                //  Now we iterate over the possible values of k2
                                //  that do not occur in the inState once k3 and k4 are removed
                                fock_t tempK2 = mask ^ (tempK3 ^ k3Occupied);
                                while(tempK2!=0)
                                {
                                    fock_t k2Occupied = tempK2 & -tempK2;      //  Isolate the right-most bit
                                    kState_t k2 = utilities::binary::HammingWeight64(k2Occupied-1);
                                    //  Get the corresponding list of k1 values
                                    iSize_t nbrK1;
                                    termTables->GetK1(kRetrieveBuffer, nbrK1, k2, k3, k4);
                                    //  Loop over the possible k1 values
                                    for(iSize_t j=0; j<nbrK1; ++j)
                                    {
                                        kState_t k1 = kRetrieveBuffer[j];
                                        fock_t k1Occupied = (fock_t)1 << k1;
                                        //  Restrict to cases that give non zero
                                        //  matrix elements, using the properties
                                        //  of spinless fermion operators
                                        //  If k1 != k3 or k4 then k1 must not occur in the 
                                        //  inState and if k2 != k3 or k4 then k2 must 
                                        //  not occur in the inState
                                        if(
                                            (!(k1Occupied & inState) || ((k1Occupied & k3Occupied) || (k1Occupied & k4Occupied))) &&
                                            (!(k2Occupied & inState) || ((k2Occupied & k3Occupied) || (k2Occupied & k4Occupied))) &&
                                            !(k1Occupied & k2Occupied)
                                        )
                                        {
                                            //  Evaluate the c^+ c^+ c c operator acting on the current state
                                            int sign;
                                            fock_t outState = utilities::Evaluate_CdCdCC(sign, k1Occupied, k2Occupied, k3Occupied, 
                                                                                         k4Occupied, inState, 
                                                                                         this->m_data.m_highestOrbital);
                                            if(0 != sign && outState>=inState)
                                            {
                                                //  If non zero then calculate the matrix element
                                                //  and find the index of the outState
                                                fock_t outIndex = this->m_fockBasis.FindFockStateIndex(mpi.m_firstTask, outState);
                                                if(outIndex < this->m_data.m_fockSpaceDim)
                                                {
                                                    T matrixElement = (double)sign*termTables->GetVkkkk(k1, k2, k3, k4);
                                                    //  Add the matrix element to the matrix
                                                    if((utilities::is_same<dcmplx,T>::value && abs(matrixElement)>_SPARSE_ZERO_TOL_) || (utilities::is_same<double, T>::value && fabs(matrixElement)>_SPARSE_ZERO_TOL_))
                                                    {
                                                        if(utilities::_DENSE_ == this->m_data.m_storageMethod)
                                                        {
                                                            this->m_fullMatrix[(i-mpi.m_firstTask)*this->m_data.m_fockSpaceDim+outIndex] += matrixElement;
                                                        }  
                                                        else if(utilities::_SPARSE_MAPPED_ == this->m_data.m_storageMethod)
                                                        {
                                                            this->m_sparseMap.Insert((i-mpi.m_firstTask), outIndex) += matrixElement;
                                                        }
                                                        else if(utilities::_SPARSE_CRS_ == this->m_data.m_storageMethod)
                                                        {
                                                            this->m_crsSparseMatrix.Insert((i-mpi.m_firstTask),outIndex)  += matrixElement;
                                                        }
                                                    }
                                                
                                                }
                                            }
                                        }
                                    }
                                    tempK2 = tempK2 ^ k2Occupied;     //  Remove right-most bit and proceed to the next iteration 
                                }
                            }
                            tempK3 = tempK3 ^ k3Occupied;     //  Remove right-most bit and proceed to the next iteration 
                        }
                        tempK4 = tempK4 ^ k4Occupied;     //  Remove right-most bit and proceed to the next iteration
                    }
                    currState = this->m_fockBasis.GetNextState(currState, it_state);
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {
                        progress.Display(i+1);
                    }
                }
                //  Finally, test to see if any matrix elements have become very small
                //  e.g. due to cancellations.
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<std::endl;
                    utilities::cout.DebuggingInfo()<<"\n\n\t- TRUNCATING MATRIX ELEMENTS WITH abs(element) < "<<_SPARSE_ZERO_TOL_<<std::endl;
                }
                if(utilities::_SPARSE_MAPPED_ == this->m_data.m_storageMethod)
                {
                    this->m_sparseMap.Trim(_SPARSE_ZERO_TOL_);
                }
                else if(utilities::_SPARSE_CRS_ == this->m_data.m_storageMethod)
                {
                    this->m_crsSparseMatrix.Trim(_SPARSE_ZERO_TOL_);
                }
                MPI_Barrier(mpi.m_comm);
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<std::endl;
                    timer.Stop();
                } 
            }
            else
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::cerr<<"\n\tERROR in Add_CdCdCC_Terms: Matrix memory not allocated!"<<std::endl;
                    mpi.m_exitFlag = true;
                }
            }
            mpi.ExitFlagTest();
	        return;
        }
    };
}   //  End namespace diagonalization
#endif
