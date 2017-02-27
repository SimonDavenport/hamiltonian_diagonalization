////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!  \file
//!		This file implements a container for a Fock basis of fermions in a 
//!     two level system
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
#include "two_level_fermion_fock_basis.hpp"

namespace diagonalization
{
    //!
    //! Comparison operator for TwoLevelFockState type
    //! used in binary search algorithm
    //!
    bool operator > (const TwoLevelFockState& state1,const TwoLevelFockState& state2)
    {
        return (state1.m_first > state2.m_first || (state1.m_first == state2.m_first && state1.m_second > state2.m_second));
    }

    //!
    //! Comparison operator for TwoLevelFockState type
    //! used in binary search algorithm
    //!
    bool operator < (const TwoLevelFockState& state1, const TwoLevelFockState& state2)
    {
        return (state1.m_first < state2.m_first || (state1.m_first == state2.m_first && state1.m_second < state2.m_second));
    }

    //!
    //! Equals operator for TwoLevelFockState type
    //! used in binary search algorithm
    //!
    bool operator == (const TwoLevelFockState& state1, const TwoLevelFockState& state2)
    {
        return (state1.m_first==state2.m_first && state1.m_second==state2.m_second);
    }

    //!
    //! Overload the out stream operator to allow simple printing of  class
    //! data.
    //!
    std::ostream& operator<< (std::ostream& os, const TwoLevelFockState& state)
    {
        os<<state.m_first<<"\t"<<state.m_second;
        return os;
    }
    
    //!
    //! Default constructor
    //!
    TwoLevelFockState::TwoLevelFockState()
    :
        m_first(0),
        m_second(0)
    {}
    
    //!
    //! Constructor for a TwoLevelFockState type
    //!
    TwoLevelFockState::TwoLevelFockState(
        const fock_t first,     //!<    First level fock state
        const fock_t second)    //!<    Second level fock state
    :
        m_first(first),
        m_second(second)
    {}

    //!
    //! Default constructor
    //!
    TwoLevelFermionFockBasis::TwoLevelFermionFockBasis()
        :
        m_basisStored(false),
        m_firstState(0)
    {}

    //!
    //! Copy constructor
    //!
    TwoLevelFermionFockBasis::TwoLevelFermionFockBasis(
        const TwoLevelFermionFockBasis& other)  //!<    A class instance to be copied
        :
        m_basis(other.m_basis),
        m_basisStored(other.m_basisStored),
        m_firstState(other.m_firstState)
    {}

    //!
    //! Copy assignment
    //!
    TwoLevelFermionFockBasis& TwoLevelFermionFockBasis::operator=(
        const TwoLevelFermionFockBasis& other)  //!<    A class instance to be copied
    {
        m_basis = other.m_basis;
        m_basisStored = other.m_basisStored;
        m_firstState = other.m_firstState;
        return *this;
    }

    //!
    //! Destructor
    //!
    TwoLevelFermionFockBasis::~TwoLevelFermionFockBasis()
    {
        this->Clear();
    }

    //!
    //! Clear current memory allocation
    //!
    void TwoLevelFermionFockBasis::Clear()
    {
        m_basis.clear();
        m_basisStored = false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to return the full two level fermion Fock space dimension
    //! 
    //! \return The dimension of the spinless fermion Fock space
    //! for the given number of particles and orbitals
    //////////////////////////////////////////////////////////////////////////////////
    fock_t TwoLevelFermionFockBasis::GetFullFockSpaceDimension(
        const iSize_t nbrParticles,         //!<    Number of particles
        const iSize_t nbrOrbitals1,         //!<    Number of first type orbitals
        const iSize_t nbrOrbitals2)         //!<    Number of second type orbitals
        const
    {
        //  Calculate the sum of the sub-Hilbert spaces obtained by partitioning
        //  the particles
        fock_t total = 0;
        for(iSize_t n0=0; n0<=nbrParticles; ++n0)
        {
            total += utilities::BinomialFromTable(nbrOrbitals1, n0)
                    *utilities::BinomialFromTable(nbrOrbitals2, nbrParticles-n0);
        }
        return total;
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to return the full two-level fermion Fock space dimension
    //! after quantum number constraints have been applied. 
    //!
    //! This is given by the cumulative dimension of the Fock space stored on
    //! the master node
    //! 
    //! \return Fock space dimension
    //////////////////////////////////////////////////////////////////////////////////
    fock_t TwoLevelFermionFockBasis::GetReducedFockSpaceDimension(
        const iSize_t nbrParticles,         //!<    Number of particles
        const iSize_t nbrOrbitals1,         //!<    Number of first type orbitals
        const iSize_t nbrOrbitals2,         //!<    Number of second type orbitals
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(m_basisStored)
        {
            fock_t fockSpaceDim = 0;
            if(mpi.m_id == 0)    //  For the master node
            {
                fockSpaceDim = m_basis.size();
            }
            mpi.Sync(&fockSpaceDim, 1, 0);
            return fockSpaceDim;
        }
        else
        {
            return this->GetFullFockSpaceDimension(nbrParticles, nbrOrbitals1, nbrOrbitals2);
        }
    }

    //!
    //! Add a single state to the basis
    //! 
    void TwoLevelFermionFockBasis::AddState(
        const fock_t state1,       //!<    First level state to add
        const fock_t state2)       //!<    Second level state to add
    {
        m_basis.push_back(TwoLevelFockState(state1, state2));
        m_basisStored = true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to determine the index in the Fock space corresponding to 
    //! a given Fock state in that space
    //!
    //! \return The index of the outState within the Fock basis. 
    //!
    //! Error code to return if the state was not found is the maximum value of
    //! the fock_t type.
    ////////////////////////////////////////////////////////////////////////////////
    void TwoLevelFermionFockBasis::FindFockStateIndex(
        const fock_t offset,            //!<    Offset in the starting index of the binary 
                                        //!     search array due e.g. to parallel data distribution
        const TwoLevelFockState& state, //!<    Two level Fock state to find index of
        fock_t& index)                  //!<    Value to be written on return for 
                                        //!     state index
        const
    {
        if(m_basisStored)
        {
            //	Determine the index of the outState using binary search
            //  Note that due to Hermiticity, we only need to binary
            //  search states above the diagonal
            //  TODO update state look up to use "Lin tables"
            index = utilities::BinarySearch<TwoLevelFockState>(state, m_basis.data(), m_basis.size());
            if(_SEARCH_ERROR_CODE_ == index)
            {
                utilities::cout.DebuggingInfo()<<"ERROR SEARCHING FOR OUT STATE!"<<std::endl;
                utilities::cout.DebuggingInfo()<<"IN STATE WAS: "<<state<<std::endl;
                return;
            }
            else
            {
                index += offset;
            }
        }
        else
        {
            std::cerr<<"ERROR in FindFockStateIndex: BASIS WAS NOT STORED "<<std::endl;
        }
    }

    //!
    //! Set the first Fock state in the list.
    //! 
    void TwoLevelFermionFockBasis::GetFirstState(
        const fock_t firstNodeIndex,            //!<    Used to check parallel
                                                //!     implementation
        std::vector<TwoLevelFockState>::const_iterator& it_state)                  
                                                //!<    Used to store
                                                //!     the current address
        const                         
    {
        if(m_basisStored)
        {
            //  Check the first node index
            if(firstNodeIndex == m_firstState)
            {
                it_state = m_basis.begin();
            }
            else
            {
                std::cerr<<"\n\tERROR in GetFirstState: non-matching first state index."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    //!
    //! Generate subsequent Fock states from the previous state
    //! 
    void TwoLevelFermionFockBasis::GetNextState(
        std::vector<TwoLevelFockState>::const_iterator& it_state)                 
                                                //!<    Used to store
                                                //!     the current address
        const
    {
        if(m_basisStored)
        {
            ++it_state;
        }
    }

    //!
    //! Get the dimension of the internal array
    //! 
    fock_t TwoLevelFermionFockBasis::GetStoredDimension() const
    {
        return m_basis.size();
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Allocate memory to store the full Fock space. Generate the Fock space 
    //! and store it in the basis data structure. 
    //!
    //! Optionally, a constraint function can be provided which will be applied
    //! to place a restriction on the allowed Fock states (for instance a quantum 
    //! number constraint)
    //!
    //! This function MUST be called on all nodes if run in parallel. When run
    //! in parallel, only the basis of states after the starting node index is
    //! stored - this is because for binary search of a Hermitian matrix element,
    //! only the upper triangular part of the matrix is required for the search.
    //!
    //! \return true if the Fock space dimension is 0, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool TwoLevelFermionFockBasis::GenerateFockSpace(
        const iSize_t nbrParticles, //!<    Number of particles
        const iSize_t nbrOrbitals1, //!<    Number of first orbitals
        const iSize_t nbrOrbitals2, //!<    Number of second orbitals 
        const std::function<bool(const TwoLevelFockState& state)>& Constraint,
                                    //!<    A constraint function to test whether the
                                    //!<    state "state" corresponds to a given quantum 
                                    //!<    number or not
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<"\n\t- BUILDING FOCK BASIS"<<std::endl;
        }
        //  Determine the set of states corresponding to the quantum number
        //  sector imposed by the constraint function.   
        //  We can call the constraint function in parallel to save time,
        //  however there will not in general be an equal distribution of
        //  allowed states over all nodes
        fock_t sectorCounter = 0;
        fock_t sectorTotalDimension = 0;
        for(iSize_t n0=0; n0<=nbrParticles; ++n0)
        {
            fock_t subspaceDim1 = utilities::BinomialFromTable(nbrOrbitals1, n0);
            fock_t subspaceDim2 = utilities::BinomialFromTable(nbrOrbitals2, nbrParticles-n0);
            mpi.DivideTasks(mpi.m_id, subspaceDim1, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
            fock_t nodeDim = mpi.m_lastTask - mpi.m_firstTask + 1; 
            //  Perform outer iteration over the leading subspace dimension
            fock_t currState1 = utilities::binary::GenerateHammingNumber(n0, mpi.m_firstTask);
            for(fock_t i=0; i<nodeDim; ++i)
            {
                //  Generate states on each node
                fock_t currState2 = utilities::binary::FirstBinaryHammingNumber64(nbrParticles-n0);
                for(fock_t i=0; i<subspaceDim2; ++i)
                {
                    if(Constraint(TwoLevelFockState(currState1, currState2)))
                    {
                        ++sectorCounter;
                    }
                    currState2 = utilities::binary::NextHammingNumber64(currState2);
                }
                currState1 = utilities::binary::NextHammingNumber64(currState1);
            }
        }
        //  Update the dimension variable to store the total number of allowed
        //  states and count the cumulative number of states on each node
        if(utilities::is_same<fock_t,uint64_t>::value)
        {
            MPI_Reduce(&sectorCounter,&sectorTotalDimension, 1, mpi.GetType<fock_t>(), MPI_SUM, 0, mpi.m_comm);
            mpi.Sync(&sectorTotalDimension, 1, 0);
        }
        else
        {
            std::cerr<<"\n\tERROR WITH MPI_Reduce on line "<<__LINE__<<" in "<<__FILE__<<": incompatible MPI data type"<<std::endl;
            exit(EXIT_FAILURE);
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<"\n\t\tFOCK BASIS CONTAINS "<<sectorTotalDimension<<" STATE(S)"<<std::endl;
        }
        //  Don't proceed if the Fock space dimension is 0 
        if(0 == sectorTotalDimension)
        {
            return true;
        }
        //  Next store a list of the allowed Fock states using the same parallelization
        fock_t* fockStateBuffer1 = 0;
        fock_t* fockStateBuffer2 = 0;
        {
            fock_t* allowedFockStates1 = 0;
            fock_t* allowedFockStates2 = 0;
            //  Allocate memory to store part of the Fock basis on each node, and
            //  on the master node allocate enough space to store the full basis
            allowedFockStates1 = new (std::nothrow) fock_t[sectorCounter];
            allowedFockStates2 = new (std::nothrow) fock_t[sectorCounter];
            fock_t* p_list1 = allowedFockStates1;
            fock_t* p_list2 = allowedFockStates2;
            for(iSize_t n0=0; n0<=nbrParticles; ++n0)
            {
                fock_t subspaceDim1 = utilities::BinomialFromTable(nbrOrbitals1, n0);
                fock_t subspaceDim2 = utilities::BinomialFromTable(nbrOrbitals2, nbrParticles-n0);
                mpi.DivideTasks(mpi.m_id, subspaceDim1, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
                fock_t nodeDim = mpi.m_lastTask - mpi.m_firstTask + 1; 
                //  Perform outer iteration over the leading subspace dimension
                fock_t currState1 = utilities::binary::GenerateHammingNumber(n0, mpi.m_firstTask);
                for(fock_t i=0; i<nodeDim; ++i)
                {
                    //  Generate states on each node
                    fock_t currState2 = utilities::binary::FirstBinaryHammingNumber64(nbrParticles-n0);
                    for(fock_t i=0; i<subspaceDim2; ++i)
                    {
                        if(Constraint(TwoLevelFockState(currState1, currState2)))
                        {
                            *p_list1 = currState1;
                            *p_list2 = currState2;
                            ++p_list1;
                            ++p_list2;
                        }
                        currState2 = utilities::binary::NextHammingNumber64(currState2);
                    }
                    currState1 = utilities::binary::NextHammingNumber64(currState1);
                }
            }
            //  MPI gather the full list on the master node
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                fockStateBuffer1 = new fock_t[sectorTotalDimension];
                fockStateBuffer2 = new fock_t[sectorTotalDimension];
            }
            MPI_Status status;
            mpi.Gather<fock_t>(allowedFockStates1, sectorCounter, fockStateBuffer1, sectorTotalDimension, 0, mpi.m_comm,status);
            mpi.Gather<fock_t>(allowedFockStates2, sectorCounter, fockStateBuffer2, sectorTotalDimension, 0, mpi.m_comm,status);
            delete[] allowedFockStates1;
            delete[] allowedFockStates2;
            //////  Sort the TwoLevelFockState in lexicographic order     //////
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                m_basis.resize(sectorTotalDimension);
                for(fock_t i=0; i<sectorTotalDimension; ++i)
                {
                    m_basis[i] = TwoLevelFockState(fockStateBuffer1[i], fockStateBuffer2[i]);
                }
                std::sort(m_basis.begin(), m_basis.end());
                //  Convert back to separate buffers for MPI send below
                for(fock_t i=0; i<sectorTotalDimension; ++i)
                {
                    fockStateBuffer1[i] = m_basis[i].m_first;
                    fockStateBuffer2[i] = m_basis[i].m_second;
                }
            }
        }
        MPI_Barrier(mpi.m_comm); 
        //  If running in parallel, then to enable the binary search of states
        //  obtained after operating with c^+c... operators, we need only store those 
        //  states that are greater than the diagonal (this is due to Hermiticity).
        //  Let's now calculate a cumulative total of the full number of rows that 
        //  will be stored on each node
        {
            mpi.DivideTasks(mpi.m_id, sectorTotalDimension, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, true);
            fock_t cumulativeDim = sectorTotalDimension - mpi.m_firstTask;
            m_firstState = mpi.m_firstTask;
            //  Allocate memory to store the allowed states
            utilities::cout.DebuggingInfo()<<"\n\t\tNODE "<<mpi.m_id<<" ALLOCATING "<<(cumulativeDim*(sizeof(fock_t))/(1024.0*1024.0));
            utilities::cout.DebuggingInfo()<<" MB TO STORE FOCK SPACE"<<std::endl;
            fock_t* allowedFockStates1 = 0;
            fock_t* allowedFockStates2 = 0;
            allowedFockStates1 = new (std::nothrow) fock_t[cumulativeDim];
            allowedFockStates2 = new (std::nothrow) fock_t[cumulativeDim];
            //  Redistribute the fockStateBuffer over the other nodes
            fock_t taskList[mpi.m_nbrProcs];
            fock_t nodeDimList[mpi.m_nbrProcs];
            MPI_Gather(&m_firstState, 1, mpi.GetType<fock_t>(), taskList, 1, mpi.GetType<fock_t>(), 0, mpi.m_comm);
            MPI_Gather(&cumulativeDim, 1, mpi.GetType<fock_t>(), nodeDimList, 1, mpi.GetType<fock_t>(), 0, mpi.m_comm);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                for(int i=1; i<mpi.m_nbrProcs; ++i)
                {
                    MPI_Send(fockStateBuffer1+taskList[i], nodeDimList[i], mpi.GetType<fock_t>(), i, 20+i, mpi.m_comm);
                    MPI_Send(fockStateBuffer2+taskList[i], nodeDimList[i], mpi.GetType<fock_t>(), i, 30+i, mpi.m_comm);
                }
                memcpy(allowedFockStates1, fockStateBuffer1, nodeDimList[0]*sizeof(fock_t));
                memcpy(allowedFockStates2, fockStateBuffer2, nodeDimList[0]*sizeof(fock_t));
                delete[] fockStateBuffer1;
                delete[] fockStateBuffer2;
            }
            else
            {
                MPI_Recv(allowedFockStates1, cumulativeDim, mpi.GetType<fock_t>(), 0, 20+mpi.m_id, mpi.m_comm, MPI_STATUS_IGNORE);
                MPI_Recv(allowedFockStates2, cumulativeDim, mpi.GetType<fock_t>(), 0, 30+mpi.m_id, mpi.m_comm, MPI_STATUS_IGNORE);
            }
            //  Convert to TwoLevelFockState list
            m_basis.resize(cumulativeDim);
            for(fock_t i=0; i<cumulativeDim; ++i)
            {
                m_basis[i] = TwoLevelFockState(allowedFockStates1[i], allowedFockStates2[i]);
            }
            delete[] allowedFockStates1;
            delete[] allowedFockStates2;
        }
        MPI_Barrier(mpi.m_comm);
        m_basisStored = true;
	    return false;
    }

    //!
    //! Return a copy of the Fock basis to the given buffer
    //!
    void TwoLevelFermionFockBasis::GetFockBasis(
        fock_t* buffer1,            //!<    Buffer to copy to
        fock_t* buffer2,            //!<    Buffer to copy to
        const fock_t dim,           //!<    Dimension of the given buffer
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
        const
    {
        if(m_basis.size()>=dim)
        {
            fock_t* p_buffer1 = buffer1;
            fock_t* p_buffer2 = buffer2;
            for(auto& it_basis : m_basis)
            {
                *buffer1 = it_basis.m_first;
                *buffer2 = it_basis.m_second;
                ++p_buffer1;
                ++p_buffer2;
            }
        }
        else
        {   
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                std::cerr<<"\n\tERROR in GetFockBasis: Given copy buffer is longer than the allocated Fock space buffer dimension"<<std::endl;
            }
            mpi.m_exitFlag = true;
        }
        mpi.ExitFlagTest();
    }

    //!
    //! Set the stored Fock basis from an external buffer
    //!
    void TwoLevelFermionFockBasis::SetFockBasis(
        fock_t* buffer1,            //!<    Buffer to copy to
        fock_t* buffer2,            //!<    Buffer to copy to
        const fock_t dim,           //!<    Dimension of the given buffer
        const fock_t firstState)    //!<    Index of first state represented
    {
        m_basis.resize(dim);
        fock_t* p_buffer1 = buffer1;
        fock_t* p_buffer2 = buffer2;
        for(fock_t i=0; i<dim; ++i, ++p_buffer1, ++p_buffer2)
        {
            m_basis[i] = TwoLevelFockState(*p_buffer1, *p_buffer2);
        }
        //  Enforce ascending Fock state order
        //  so that we can populate only the upper triangular part of
        //  a symmetric/Hermitian matrix
        std::sort(m_basis.begin(), m_basis.end());
        m_firstState = firstState;
        m_basisStored = true;
    }
}   //  End namespace diagonalization
