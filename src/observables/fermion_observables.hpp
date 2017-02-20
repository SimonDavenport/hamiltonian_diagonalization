////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 22/03/2015
//!
//!  \file
//!		This file contains a list of functions to perform analysis on 
//!     eigenvectors obtained from a diagonalized Hamiltonian
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

#ifndef _OBSERVABLES_HPP_INCLUDED_
#define _OBSERVABLES_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../basis/fermion_fock_basis.hpp"
#include "../utilities/mathematics/fermions.hpp"          
                                                //  For binary encoding of fermion
                                                //  operator commutations
#include "../utilities/algorithms/quick_sort.hpp" 
                                                //  For quick sort algorithm

#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{

struct Observables
{
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the total probability associated with each 
    //! orbital for a given eigenvector i.e. a table of < c^+_k c_k >
    //!
    //! This function MUST be called on all nodes if run in parallel
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    template <class H>
    void DensityFunction(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to analyse
        double* probabilityList,                //!<    Pre-allocated array of size m_nbrOrbitals
                                                //!     To store a map of occupation probabilities    
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in DensityFunction: eigenvectorIndex out of range."<<std::endl;
                }
                
                return;
            }

            //  Ensure that the probability list is initialized to zero
            
            for(iSize_t i = 0;i<hamiltonian.m_data.m_nbrOrbitals;++i)
            {
                probabilityList[i] = 0.0;
            }

            //  Loop over the Fock basis. For each state, pick out the N occupied orbitals
            //  and update the corresponding bin in the probability table with the
            //  value stored in the chosen eigenvector

            mpi.DivideTasks(mpi.m_id,hamiltonian.m_data.m_fockSpaceDim,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);

            std::vector<fock_t>::const_iterator it_state;
            fock_t currState = hamiltonian.m_fockBasis.GetFirstState(mpi.m_firstTask,hamiltonian.m_data.m_nbrParticles,it_state);

            auto it_eigs = hamiltonian.m_eigenvectors.begin() + eigenvectorIndex*hamiltonian.GetEigenvectorNodeDim();

            for(iSize_t i = mpi.m_firstTask;i<=mpi.m_lastTask;++i,++it_eigs)
            {
                fock_t inState = currState;
            
                double currProb = std::real(*it_eigs * std::conj(*it_eigs));

                //std::cout<<"state: "<<(std::bitset<64>(inState))<<" prob: "<<currProb<<std::endl;

                while(inState!=0)
                {
                    //  Isolate successive right-most bits until we've included all
                    //  occupied states in the current basis state
                    
                    //  Isolate the right-most bit
                    fock_t kOccupied = inState & -inState;      

                    //  Bin the probability
                    probabilityList[utilities::binary::HammingWeight64(kOccupied-1)] += currProb;
                    
                    //  Remove right-most bit and proceed to the next iteration
                    inState = inState ^ kOccupied;     
                }
                
                currState = hamiltonian.m_fockBasis.GetNextState(currState,it_state);
            }

            //  Sum up all the contributions to the probability list
            //  and keep the total on the master node
            
            double* buffer = new double[hamiltonian.m_data.m_nbrOrbitals];
            
            MPI_Status status;
            
            mpi.Reduce<double>(probabilityList,buffer,hamiltonian.m_data.m_nbrOrbitals,0,mpi.m_comm,status);
            
            delete[] buffer;

            MPI_Barrier(mpi.m_comm);
        }
        else
        {
            exit(EXIT_FAILURE);
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the total density-density operator associated with 
    //! a pair of orbitals for a given eigenvector i.e. a table of 
    //! < c^+_k c_k c^+ q c_q>
    //!
    //! This function MUST be called on all nodes if run in parallel
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <class H>
    void DensityDensityFunction(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to analyse
        double* probabilityList,                //!<    Pre-allocated array of size
                                                //!     m_nbrOrbitals*m_nbrOrbitals
                                                //!     to store a map of occupation probabilities
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in DensityDensityFunction: eigenvectorIndex out of range."<<std::endl;
                }
                
                return;
            }

            //  Ensure that the probability list is initialized to zero
            
            for(iSize_t i = 0;i<hamiltonian.m_data.m_nbrOrbitals*hamiltonian.m_data.m_nbrOrbitals;++i)
            {
                probabilityList[i] = 0.0;
            }

            //  Apply the operator c^+_k c_k c^+ q c_q to each Fock state
            //  and then determine the matrix element with the initial
            //  Fock state
            
            //  A non-zero result will occur in only two cases:
            //  
            //  1. The initial state contains q, but not k. This requires k=q
            //  (and gives the value n_k^2)
            //
            //  2. The initial state contains both q and k. This requires k!=q

            //  Loop over the Fock basis. For each state, pick out the N occupied orbitals
            //  corresponding to q. We use this to calculate n_k^2, then we loop over the 
            //  remaining N-1 orbitals containing the possible values of k. The value
            //  of n_k n_q is added to a bin (k,q)

            mpi.DivideTasks(mpi.m_id,hamiltonian.m_data.m_fockSpaceDim,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);

            std::vector<fock_t>::const_iterator it_state;
            fock_t currState = hamiltonian.m_fockBasis.GetFirstState(mpi.m_firstTask,hamiltonian.m_data.m_nbrParticles,it_state);

            auto it_eigs = hamiltonian.m_eigenvectors.begin() + eigenvectorIndex*hamiltonian.GetEigenvectorNodeDim();

            for(iSize_t i = mpi.m_firstTask;i<=mpi.m_lastTask;++i,++it_eigs)
            {
                fock_t inState = currState;
            
                double currProb = std::real(*it_eigs * std::conj(*it_eigs));

                //std::cout<<"state: "<<(std::bitset<64>(inState))<<" prob: "<<currProb<<std::endl;

                fock_t tempState = inState;

                while(tempState!=0) //  Iterate over N terms
                {
                    //  Isolate successive right-most bits until we've included all
                    //  occupied states in the current basis state
                    
                    //  Isolate the right-most bit
                    fock_t kOccupied = tempState & -tempState;      

                    //  Bin the probability as n_k n_k
                    
                    uint64_t index = utilities::binary::HammingWeight64(kOccupied-1);
                    
                    probabilityList[index*hamiltonian.m_data.m_nbrOrbitals+index] += currProb;
                    
                    //  Remove right-most bit and proceed to the next iteration
                    tempState = tempState ^ kOccupied;
                    
                    //  Remove only the single selected bit from the in state,
                    //  leaving N-1 set bits 
                    
                    fock_t oneRemoved = inState ^ kOccupied;  
                    
                    //  Inner iteration over all but the removed bit
                    
                    while(oneRemoved!=0)    //  Iterate over N-1 terms
                    {
                        //  Isolate the right-most bit
                        fock_t qOccupied = oneRemoved & -oneRemoved;
                        
                        //  Bin the probability as n_k n_q
                        
                        probabilityList[index*hamiltonian.m_data.m_nbrOrbitals+utilities::binary::HammingWeight64(qOccupied-1)] += currProb;
                    
                        //  Remove right-most bit and proceed to the next iteration
                        oneRemoved = oneRemoved ^ qOccupied;
                    }
                }
                
                currState = hamiltonian.m_fockBasis.GetNextState(currState,it_state);
            }

            MPI_Barrier(mpi.m_comm);
            
            //  Allocate a temporary buffer for MPI communication
            
            double* reduceBuffer = new double[hamiltonian.m_data.m_nbrOrbitals*hamiltonian.m_data.m_nbrOrbitals];
            
            MPI_Status status;
            
            mpi.Reduce(probabilityList,reduceBuffer,hamiltonian.m_data.m_nbrOrbitals*hamiltonian.m_data.m_nbrOrbitals,0,mpi.m_comm,status);
            
            //  MPI_sync the return array over all nodes
            
            mpi.Sync(probabilityList,hamiltonian.m_data.m_nbrOrbitals*hamiltonian.m_data.m_nbrOrbitals,0);

            delete[] reduceBuffer;
        }
        else
        {
            exit(EXIT_FAILURE);
        }
        
        return;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the total translational density-density operator 
    //! associated with  a pair of orbitals for a given eigenvector
    //!
    //! Given by:
    //!
    //! S(q) = 1/N sum_{k1,k2} <c^+_k2-q c_k2 c^+_k1+q c_k1>
    //!
    //! Where q takes any value in the simulation grid
    //!
    //! This function MUST be called on all nodes if run in parallel
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    template <class H>
    void TranslationalDensityDensityFunction(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to analyse
        std::function<kState_t(const kState_t k1,const kState_t k2)> StateAddition,
                                                //!<    A function to calculate the output from
                                                //!     summing any two k-states k1+k2
        std::function<kState_t(const kState_t k1,const kState_t k2)> StateSubtraction,          
                                                //!<    A function to calculate the output from
                                                //!     subtracting any two k-states k1-k2
        std::vector<double>& functionList,      //!<    Pre-allocated array to store a map of 
                                                //!     occupation probabilities for each G value
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in TranslationalDensityDensityFunction: eigenvectorIndex out of range."<<std::endl;
                }
                
                return;
            }

            //  Ensure that the function list is initialized to zero

            for(auto& it : functionList)
            {
                it = 0.0;
            }

            //  Apply the operator c^+_k2-G c_k2 c^+_k1+G c_k1 to each Fock state
            //  and then determine the matrix element with the initial Fock state

            //  Loop over the Fock basis

            mpi.DivideTasks(mpi.m_id,hamiltonian.m_data.m_fockSpaceDim,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);

            std::vector<fock_t>::const_iterator it_state;
            fock_t currState = hamiltonian.m_fockBasis.GetFirstState(mpi.m_firstTask,hamiltonian.m_data.m_nbrParticles,it_state);

            auto it_eigs = hamiltonian.m_eigenvectors.begin() + eigenvectorIndex*hamiltonian.GetEigenvectorNodeDim();
            
            for(iSize_t i = mpi.m_firstTask;i<=mpi.m_lastTask;++i,++it_eigs)
            {
                fock_t inState = currState;
            
                double currProb = std::real(*it_eigs * std::conj(*it_eigs));

                //std::cout<<"state: "<<(std::bitset<64>(inState))<<" prob: "<<currProb<<std::endl;

                fock_t tempState = inState;

                while(tempState!=0) //  Iterate over N terms
                {
                    //  Isolate successive right-most bits until we've included all
                    //  occupied states in the current basis state
                    
                    //  Isolate the right-most bit
                    fock_t k1Occupied = tempState & -tempState;
                    
                    kState_t k1 = utilities::binary::HammingWeight64(k1Occupied-1);

                    //PRINT("k1",k1);
                    
                    for(kState_t q = 0;q<functionList.size();++q)
                    {
                        //PRINT("q",q);
                    
                        fock_t k1PlusQOccupied = (fock_t)1 << StateAddition(k1,q);
                        
                        //PRINT("k1+q",StateAddition(q ,k1));
                        
                        //PRINT("k1PlusQOccupied",(std::bitset<64>(k1PlusQOccupied)));
                        
                        //  Act on the current in state with CdC
                        
                        int sign;

                        fock_t outState = utilities::Evaluate_CdC(sign,k1PlusQOccupied,k1Occupied,inState,hamiltonian.m_data.m_highestOrbital);
                        
                        //PRINT("outState",(std::bitset<64>(outState)));
                        
                        if(sign!=0)
                        {
                            fock_t tempState2 = outState;
                            
                            while(tempState2!=0) //  Iterate over N terms
                            {
                                //  Isolate the right-most bit
                                fock_t k2Occupied = tempState2 & -tempState2;
                            
                                kState_t k2 = utilities::binary::HammingWeight64(k2Occupied-1);
                            
                                //PRINT("k2",k2);
                            
                                //PRINT("k2Occupied",(std::bitset<64>(k2Occupied)));
                            
                                //  Isolate successive right-most bits until we've included all
                                //  occupied states in the current basis state

                                //PRINT("k2-q",StateSubtraction(k2,q));

                                fock_t k2MinusQOccupied = (fock_t)1 << StateSubtraction(k2,q);
                                
                                //PRINT("k2MinusQOccupied",(std::bitset<64>(k2MinusQOccupied)));
                                
                                //  Act on the current in state with CdC
                                
                                int sign2;

                                fock_t outState2 = utilities::Evaluate_CdC(sign2,k2MinusQOccupied,k2Occupied,outState,hamiltonian.m_data.m_highestOrbital);
                                
                                //PRINT("outState2",(std::bitset<64>(outState2)));
                                
                                //PAR_PAUSE();
                                
                                //  Include the result if the out state matches the in state
                                
                                if(outState2 == inState && sign2 != 0)
                                {
                                    functionList[q] += sign*sign2*currProb;
                                }
                                
                                //  Remove right-most bit and proceed to the next iteration
                                tempState2 = tempState2 ^ k2Occupied;
                            }
                        }
                    }
                    
                    //  Remove right-most bit and proceed to the next iteration
                    tempState = tempState ^ k1Occupied;
                }
                
                currState = hamiltonian.m_fockBasis.GetNextState(currState,it_state);
            }
            
            MPI_Barrier(mpi.m_comm);

            //  Allocate a temporary buffer for MPI communication
            
            double* reduceBuffer = new double[functionList.size()];
            
            MPI_Status status;
            
            mpi.Reduce(functionList.data(),reduceBuffer,functionList.size(),0,mpi.m_comm,status);
            
            //  MPI_sync the return array over all nodes
            
            mpi.Sync(functionList.data(),functionList.size(),0);

            delete[] reduceBuffer;
            
            //  Divide all values by a factor of N
            
            for(auto& it : functionList)
            {
                it /= hamiltonian.m_data.m_nbrParticles;
            }
        }
        else
        {
            exit(EXIT_FAILURE);
        }
        
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the total rotational density-density operator 
    //! associated with  a pair of orbitals for a given eigenvector
    //!
    //! Given by:
    //!
    //! R(m) = 1/6N sum_{n,k} exp(i*m*n*pi/3) <c^+_{R_n k} c_{R_n k} c^+_k c_k>
    //!
    //! Where m and n take values in {0,1,2,3,4,5}
    //!
    //! This function MUST be called on all nodes if run in parallel
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    template <class H>
    void RotationalDensityDensityFunction(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to analyse
        std::function<kState_t(const kState_t k,const int nbrApplied)> R60,
                                                //!<    A function to apply the R60
                                                //!     rotation on the given orbital
        std::function<kState_t(const kState_t k,const int nbrApplied)> InverseR60,
                                                //!<    A function to apply the inverse R60
                                                //!     rotation on the given orbital
        dcmplx* functionList,                   //!<    Pre-allocated array of size 6
                                                //!     to store a map of occupation probabilities
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        const double tol = 0.00000000001;
        const int nRot   = 6;
    
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in RotationalDensityDensityFunction: eigenvectorIndex out of range."<<std::endl;
                }
                
                return;
            }

            //  Ensure that the function list is initialized to zero
            
            for(int delta = 0;delta<nRot;++delta)
            {
                functionList[delta] = 0.0;
            }

            //  For each Fock state, determine the set of k1,...,k4 and alpha,alpha'
            //  that will return the original state
            
            //  Loop over the Fock basis

            mpi.DivideTasks(mpi.m_id,hamiltonian.m_data.m_fockSpaceDim,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);

            std::vector<fock_t>::const_iterator it_state;
            fock_t currState = hamiltonian.m_fockBasis.GetFirstState(mpi.m_firstTask,hamiltonian.m_data.m_nbrParticles,it_state);

            auto it_eigs = hamiltonian.m_eigenvectors.begin() + eigenvectorIndex*hamiltonian.GetEigenvectorNodeDim();
            
            //fock_t mask = hamiltonian.m_data.m_highestOrbital - 1;
            
            for(iSize_t i = mpi.m_firstTask;i<=mpi.m_lastTask;++i,++it_eigs)
            {
                fock_t inState = currState;
            
                double currProb = std::real(*it_eigs * std::conj(*it_eigs));

                //std::cout<<"state: "<<(std::bitset<64>(inState))<<" prob: "<<currProb<<std::endl;
                
                if(fabs(currProb) > tol)
                {
                    fock_t tempState1 = inState;
                    
                    while(tempState1!=0) //  Iterate over N terms
                    {
                        //  Isolate successive right-most bits until we've included all
                        //  occupied states in the current basis state
                        
                        //  Isolate the right-most bit
                        fock_t k1Occupied = tempState1 & -tempState1;
                        
                        kState_t k1 = utilities::binary::HammingWeight64(k1Occupied-1);
       
                        fock_t tempState2 = inState;
                        
                        while(tempState2!=0) //  Iterate over N terms
                        {
                            //  Isolate the right-most bit
                            fock_t k2TempOccupied = tempState2 & -tempState2;
                            
                            kState_t k2Temp = utilities::binary::HammingWeight64(k2TempOccupied-1);
                            
                            //  Sum over alpha that will produce this k2Temp value from k1
                            
                            for(int alpha = 0;alpha<nRot;++alpha)
                            {
                                kState_t k2 = InverseR60(k2Temp,alpha);
                                    
                                if(k2 == k1)
                                {
                                    for(int m = 0;m<nRot;++m)
                                    {
                                        functionList[m] += currProb*exp(I*(PI/3.0)*(double)(m*alpha));
                                    }
                                }
                            }
                                
                            //  Remove right-most bit and proceed to the next iteration
                            tempState2 = tempState2 ^ k2TempOccupied;
                        }

                        //  Remove right-most bit and proceed to the next iteration
                        tempState1 = tempState1 ^ k1Occupied;
                    }
                }
                
                currState = hamiltonian.m_fockBasis.GetNextState(currState,it_state);
            }

            MPI_Barrier(mpi.m_comm);
            
            //  Allocate a temporary buffer for MPI communication
            
            dcmplx reduceBuffer[nRot];
            
            MPI_Status status;
            
            mpi.Reduce(functionList,reduceBuffer,nRot,0,mpi.m_comm,status);
            
            //  MPI_sync the return array over all nodes
            
            mpi.Sync(functionList,nRot,0);
            
            //PRINT("functionList",functionList,nRot);
            
            //  Divide all values by a factor of nRot*N
            
            for(int delta = 0;delta<nRot;++delta)
            {
                functionList[delta] /= nRot*hamiltonian.m_data.m_nbrParticles;
            }
        }
        else
        {
            exit(EXIT_FAILURE);
        }
        
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function returns the nbrStates most probable eigenstates and
    //! their amplitudes to a set of buffers.
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <class H>
    void GetMostProbable(
        const H& hamiltonian,               //!<    Hamiltonian object
        const iSize_t nbrStates,            //!<    The number of most probable states
                                            //!     to be returned
        fock_t* states,                     //!<    A buffer to store the returned Fock
                                            //!     states ( of size nbrStates*nbr eigenvectors)
        dcmplx* amplitudes,                 //!<    A buffer to store the associated amplitudes
                                            //!     (of size nbrStates*nbr eigenvectors)
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {
            //  Note that we cannot return more most probable amplitudes than the fock space dimension

            iSize_t actualNbrStates = std::min(hamiltonian.m_data.m_fockSpaceDim,(fock_t)nbrStates);
            
            //  Note that we cannot sort more amplitudes than are stored on a given node
            
            iSize_t nbrStatesPerNode = std::min(hamiltonian.GetEigenvectorNodeDim(),(fock_t)nbrStates);
            
            iSize_t totalNbrStates = 0;
            
            //  Get the total number of states sorted on each node
                
            MPI_Allreduce(&nbrStatesPerNode,&totalNbrStates,1,mpi.GetType<iSize_t>(),MPI_SUM,mpi.m_comm);
            
            for(iSize_t e=0;e<hamiltonian.m_data.m_nbrEigenvalues;++e)
            {
                if(0 == mpi.m_id)    //  FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\n\t For eigenvector "<<e<<std::endl;
                }

                //  On each node make a copy of the eigenvector data and the Fock
                //  state data into buffers to be used for sorting
                
                dcmplx* eigsBuffer = 0;
                fock_t* fockBuffer = 0;
                fock_t eigenvectorNodeDim = hamiltonian.GetEigenvectorNodeDim();
 
                if(hamiltonian.GetEigenvectorNodeDim()>0)
                {
                    eigsBuffer = new dcmplx[eigenvectorNodeDim];
                    fockBuffer = new fock_t[eigenvectorNodeDim];

                    memcpy(eigsBuffer,&hamiltonian.m_eigenvectors[0]+e*eigenvectorNodeDim,eigenvectorNodeDim*sizeof(dcmplx));
                    
                    hamiltonian.m_fockBasis.GetFockBasis(fockBuffer,eigenvectorNodeDim,mpi);
                    
                    //  Partially sort each eigenvector in place into descending order 
                    //  (using the abs value as a comparison function)
                    
                    utilities::PartialQuickSort<dcmplx,fock_t,_DESCENDING_ORDER_>(eigsBuffer,fockBuffer,eigenvectorNodeDim,nbrStatesPerNode);
                }

                //  Gather the sorted values onto the master node, then perform a further
                //  sort to single out the overall highest nbrStates values to be returned
                
                //  Note that we only need to use the first nbrStates values from each node 
                
                dcmplx* gatherAmplitudesBuffer = 0;
                fock_t* gatherStatesBuffer     = 0;
                
                if(0 == mpi.m_id)    //  FOR THE MASTER NODE
                {
                    gatherAmplitudesBuffer = new dcmplx[nbrStates*mpi.m_nbrProcs];
                    gatherStatesBuffer     = new fock_t[nbrStates*mpi.m_nbrProcs];
                }

                MPI_Status status;

                mpi.Gather(eigsBuffer,nbrStatesPerNode,gatherAmplitudesBuffer,
                           totalNbrStates,0,mpi.m_comm,status);
                
                mpi.Gather(fockBuffer,nbrStatesPerNode,gatherStatesBuffer,
                           totalNbrStates,0,mpi.m_comm,status);

                //  Now we can deallocate the eigenvectors/fock states buffers on each node

                delete[] eigsBuffer;
                delete[] fockBuffer;

                //  Partially sort the combined list of values to obtain the most probable
                //  states, then copy these into the return vectors
                
                if(0 == mpi.m_id)    //  FOR THE MASTER NODE 
                {
                    utilities::PartialQuickSort<dcmplx,fock_t,_DESCENDING_ORDER_>(gatherAmplitudesBuffer,gatherStatesBuffer,totalNbrStates,actualNbrStates);

                    memcpy(amplitudes+e*actualNbrStates,gatherAmplitudesBuffer,nbrStates*sizeof(dcmplx));
                    memcpy(states+e*actualNbrStates,gatherStatesBuffer,nbrStates*sizeof(fock_t));
                    
                    for(iSize_t i=0;i<nbrStates;++i)
                    {
                        utilities::cout.AdditionalInfo()<<"\t State "<<std::setw(15)<<std::left<<states[i+e*actualNbrStates]<<"\t Amplitude "<<std::setw(15)<<std::left<<amplitudes[i+e*actualNbrStates]<<std::endl;
                    }
                }

                //  Free up the gather buffer memory

                if(0 != gatherAmplitudesBuffer) delete[] gatherAmplitudesBuffer;          
                if(0 != gatherStatesBuffer)     delete[] gatherStatesBuffer;
            }
        }
        
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the fractional participation ratio and store it in a file
    //!
    //! Fractional participation ratio is defined as 1/ D sum |psi|^4.
    //! This function gives a measure of the proportion of orbitals over
    //! which the many body wave function may be considered to be
    //! localized
    //!
    //! \return Participation ratio
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <class H>
    double GetParticipationRatio(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to analyse
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in GetParticipationRatio: eigenvectorIndex out of range."<<std::endl;
                }
                
                return 0.0;
            }
        
            //  Divide up calculation over multiple nodes

            mpi.DivideTasks(mpi.m_id,hamiltonian.m_data.m_fockSpaceDim,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);
            
            auto it_eigs = hamiltonian.m_eigenvectors.begin() + eigenvectorIndex*hamiltonian.GetEigenvectorNodeDim();

            double participationRatio = 0.0;
            double totalParticipationRatio = 0.0;

            for(iSize_t i = mpi.m_firstTask;i<=mpi.m_lastTask;++i,++it_eigs)
            {
                participationRatio += std::pow(abs(*it_eigs),4.0);
            }
            
            MPI_Barrier(mpi.m_comm);
            
            //  Sum the participation ratio over all nodes
            
            MPI_Allreduce(&participationRatio,&totalParticipationRatio,1,mpi.GetType<double>(),MPI_SUM,mpi.m_comm);
            
            //  An estimate of the proportion of Fock states that
            //  mostly contribute to the eigenvector is given by
            //  1/ Dsum |psi|^4:

            return 1.0/(totalParticipationRatio);   //hamiltonian.m_data.m_fockSpaceDim
        }
        else
        {
            exit(EXIT_FAILURE);
        }
        
        return 0.0;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the particle entanglement spectrum and store in a file
    //!
    //! TODO implement particle entanglement spectrum
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <class H>
    void ParticleEntanglementSpectrumToFile(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to decompose
        const iSize_t nbrA,                     //!<    Number of particles in the A cut
        const std::string fileName,             //!<    Output file name
        const utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        {   
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in ParticleEntanglementSpectrumToFile: eigenvectorIndex out of range."<<std::endl;
                }
                
                return;
            }

            //  Calculate the particle entanglement spectrum
            
            //  Write the results to a text file

            if(0 == mpi.m_id)
            {
            
            
            }
            
            exit(EXIT_FAILURE);
        }
        else
        {
            exit(EXIT_FAILURE);
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the orbital entanglement spectrum and store in a file
    //!
    //! The OES is calculate directly from the eigenvector representing our
    //! wave function. First we group the set of site indices according to
    //! our particular choice of orbital cut. Then we take the SVD of the 
    //! resulting (dense) matrix.
    //!
    //! TODO implement orbital entanglement spectrum
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    template <class H>
    void OrbitalEntanglementSpectrumToFile(
        const H& hamiltonian,                   //!<    Hamiltonian object
        const iSize_t eigenvectorIndex,         //!<    Index of eigenvector to decompose
        const iSize_t orbitalCut,               //!<    Orbital label at the cut
        const std::string fileName,             //!<    Output file name#
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
    /*
        if(hamiltonian.m_data.m_eigenvectorsCalculated)
        { 
            if(eigenvectorIndex>hamiltonian.m_data.m_nbrEigenvalues)
            {
                if(0 == mpi.m_id)
                {
                    std::cerr<<"\n\tERROR in OrbitalEntanglementSpectrumToFile: eigenvectorIndex out of range."<<std::endl;
                }
                
                return;
            }

            //////      Gather the distributed eigenvector      ////////////////////
            
            dcmplx* gatheredEigenvector = 0;
            
            if(0 == mpi.m_id)    //  For the master node
            {   
                //  Allocate space for the reduced density matrix on the master node
                gatheredEigenvectorBuffer = new dcmplx[hamiltonian.m_data.m_fockSpaceDim];
            }

            MPI_Status status;

            fock_t eigenvectorNodeDim = hamiltonian.GetEigenvectorNodeDim();
            mpi.Gather(&m_eigenvectors[0]+eigenvectorIndex*eigenvectorNodeDim,eigenvectorNodeDim,gatheredEigenvectorBuffer,hamiltonian.m_data.m_fockSpaceDim,0,mpi.m_comm,status);

            if(0 == mpi.m_id)    //  For the master node
            {          
                //////      Initialize parameters for the OES calculation       ////
            
                utilities::IndexManager manager();  //  Index manager class to support
                                                //  arbitrary index rearrangements 
                
                double* entanglementEigs;       //  Buffer to store the singular values
                fock_t rank;           //  Rank of the reduced density matrix

                manager.SetSystemData(hamiltonian.hamiltonian.m_data.m_nbrParticles,hamiltonian.m_data.m_nbrOrbitals);
                manager.SetVectorData(gatheredEigenvectorBuffer,hamiltonian.m_data.m_fockSpaceDim);
                manager.SetOrbitalsData(hamiltonian.m_fockBasis.GetStartStatePtr(),hamiltonian.m_data.m_fockSpaceDim);

                //  Release gathered eigenvector buffer
                delete[] gatheredEigenvectorBuffer;

                //  Perform an in-place reordering of the reduced density matrix
                //  indices specified by the orbital cut
                
                manager.Rearrange(orbitalCut);
                
                //  Perform a singular value decomposition of the reduced density
                //  matrix in order to get the entanglement eigenvectors

                manager.GetRank(&rank);

                entanglementEigs = new double[rank];

                manager.SingularValueDecomposition(entanglementEigs);

                //  Write the results to a text file

                std::ofstream f_out;
                
                f_out.open(fileName.str().c_str(), std::ios::out);
                
                if(!f_out.is_open())
                {
                    std::cerr<<"\n\tERROR in OrbitalEntanglementSpectrumToFile: Cannot open file "<<fileName.str()<<std::endl;
                    
                    mpi.m_exitFlag = true;
                }
                else
                {
                    f_out.precision(15);
                
                    for(fock_t i=0;i<rank;i++)
                    {
                        f_out<<entanglementEigs[i]<<"\n";
                    }

                    f_out.close();
                    
                    utilities::cout.SecondaryOutput()<<"\n\tOrbital entanglement spectrum hamiltonian.m_data output to file "<<fileName.str()<<std::endl;
                }
                
                delete[] entanglementEigs;
            }
            
            mpi.ExitFlagTest();
        }
        else
        {
            exit(EXIT_FAILURE);
        }
        
        return;
        
    */
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
};

}   //  End namespace diagonalization 

#endif
