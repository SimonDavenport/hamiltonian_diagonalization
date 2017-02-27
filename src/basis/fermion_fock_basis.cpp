////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file implements a container for a Fock basis of fermions
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
#include "fermion_fock_basis.hpp"

namespace diagonalization
{
    //!
    //! Default constructor
    //!
    FermionFockBasis::FermionFockBasis()
        :
        m_basisStored(false),
        m_firstState(0)
    {}

    //!
    //! Copy constructor
    //!
    FermionFockBasis::FermionFockBasis(
        const FermionFockBasis& other)  //!<    A class instance to be copied
        :
        m_basis(other.m_basis),
        m_basisStored(other.m_basisStored),
        m_firstState(other.m_firstState)
    {
    }

    //!
    //! Copy assignment
    //!
    FermionFockBasis& FermionFockBasis::operator=(
        const FermionFockBasis& other)  //!<    A class instance to be copied
    {
        m_basis = other.m_basis;
        m_basisStored = other.m_basisStored;
        m_firstState = other.m_firstState;
        return *this;
    }

    //!
    //! Destructor
    //!
    FermionFockBasis::~FermionFockBasis()
    {
        this->Clear();
    }

    //!
    //! Clear current memory allocation
    //!
    void FermionFockBasis::Clear()
    {
        m_basis.clear();
        m_basisStored = false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to return the full spinless fermion Fock space dimension
    //! 
    //! \return The dimension of the spinless fermion Fock space
    //! for the given number of particles and orbitals
    //////////////////////////////////////////////////////////////////////////////////
    fock_t FermionFockBasis::GetFullFockSpaceDimension(
        const iSize_t nbrParticles,     //!<    Number of particles
        const iSize_t nbrOrbitals)      //!<    Number of orbitals
        const
    {
        return utilities::BinomialFromTable(nbrOrbitals, nbrParticles);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to return the spinless fermion Fock space dimension
    //! after quantum number constraints have been applied. 
    //!
    //! This is given by the cumulative dimension of the Fock space stored on
    //! the master node
    //! 
    //! \return Fock space dimension
    //////////////////////////////////////////////////////////////////////////////////
    fock_t FermionFockBasis::GetReducedFockSpaceDimension(
        const iSize_t nbrParticles,         //!<    Number of particles
        const iSize_t nbrOrbitals,          //!<    Number of orbitals
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        fock_t fockSpaceDim = 0;
        if(mpi.m_id == 0)    //  For the master node
        {
            fockSpaceDim = m_basis.size();
        }
        mpi.Sync(&fockSpaceDim, 1, 0);
        if(m_basisStored)
        {
            return fockSpaceDim;
        }
        else
        {
            return this->GetFullFockSpaceDimension(nbrParticles, nbrOrbitals);
        }
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
    fock_t FermionFockBasis::FindFockStateIndex(
        const fock_t offset,       //!<    Offset in the starting index of the binary 
                                   //!     search array due e.g. to parallel data distribution
        const fock_t state)        //!<    Fock state to find the index of
        const
    {
        if(m_basisStored)
        {
            //	Determine the index of the outState using binary search. Note that 
            //  due to Hermiticity, we only need to binary search states above the diagonal
            const fock_t tempIndex = utilities::BinarySearch<fock_t>(state, &m_basis[0], m_basis.size());
            if(_SEARCH_ERROR_CODE_ == tempIndex)
            {
                utilities::cout.DebuggingInfo()<<"ERROR SEARCHING FOR OUT STATE!"<<std::endl;
                utilities::cout.DebuggingInfo()<<"STATE WAS: "<<state<<std::endl;
                return tempIndex;
            }
            else
            {
                return offset + tempIndex;
            }
        }
        else
        {
            std::cerr<<"ERROR in FindFockStateIndex: BASIS WAS NOT STORED "<<std::endl;
            return _SEARCH_ERROR_CODE_;
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief Set the first Fock state in the list. If the list is not 
    //! populated, use a combinatorial generation scheme, which applies
    //! when there are no quantum number constraints
    //////////////////////////////////////////////////////////////////////////////////
    fock_t FermionFockBasis::GetFirstState(
        const fock_t firstNodeIndex,            //!<    Optionally used in
                                                //!     combinatorial scheme
        const iSize_t nbrParticles,             //!<    Optionally used in 
                                                //!     combinatorial scheme
        std::vector<fock_t>::const_iterator& it_state)                   
                                                //!<    Optionally used to store
                                                //!     the current address
        const                                        
    {
        if(m_basisStored)
        {
            if(firstNodeIndex == m_firstState)
            {
                it_state = m_basis.begin();
                return m_basis[0];
            }
            else
            {
                std::cerr<<"\n\tERROR in GetFirstState: non-matching first state index."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            return utilities::binary::GenerateHammingNumber(nbrParticles,firstNodeIndex);
        } 
    }

    //!
    //! Generate subsequent Fock states from the previous state
    //! 
    fock_t FermionFockBasis::GetNextState(
        const fock_t state,                     //!<    Current state value
        std::vector<fock_t>::const_iterator& it_state)                   
                                                //!<    Optionally used to store
                                                //!     the current address
        const
    {
        if(m_basisStored)
        {
            ++it_state;
            return *(it_state);
        }
        else
        {
            return utilities::binary::NextHammingNumber64(state);
        }
    }

    //!
    //! Get the dimension of the internal array
    //! 
    fock_t FermionFockBasis::GetStoredDimension() const
    {
        return m_basis.size();
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Allocate memory to store the full Fock space. Generate the Fock space 
    //! and store it in the m_basis data structure. 
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
    bool FermionFockBasis::GenerateFockSpace(
        const iSize_t nbrParticles, //!<    Number of particles
        const iSize_t nbrOrbitals,  //!<    Number of orbitals 
        const std::function<bool(const fock_t state)>& Constraint,
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
        fock_t fullFockSpaceDim = this->GetFullFockSpaceDimension(nbrParticles,nbrOrbitals);
        mpi.DivideTasks(mpi.m_id, fullFockSpaceDim, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
        fock_t fullNodeDim = mpi.m_lastTask - mpi.m_firstTask + 1;
        fock_t sectorCounter = 0;
        fock_t sectorTotalDimension = 0;
        {
            //  Generate first state on each node
            fock_t currState = utilities::binary::GenerateHammingNumber(nbrParticles,mpi.m_firstTask);
            //  Count the set of states in quantum number sector specified by the Constraint function
            for(fock_t i=0; i<fullNodeDim; ++i)
            {
                if(Constraint(currState))
                {
                    ++sectorCounter;
                }
                currState = utilities::binary::NextHammingNumber64(currState);
            }
            //  Update the dimension variable to store the total number of allowed
            //  states and count the cumulative number of states on each node
            MPI_Reduce(&sectorCounter, &sectorTotalDimension, 1, mpi.GetType<fock_t>(), MPI_SUM, 0, mpi.m_comm);
            mpi.Sync(&sectorTotalDimension, 1, 0);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t\tFOCK BASIS CONTAINS "<<sectorTotalDimension<<" STATE(S)"<<std::endl;
            }
            //  Don't proceed if the Fock space dimension is 0
            if(0 == sectorTotalDimension)
            {
                return true;
            }
        }
        //  Next store a list of the allowed Fock states using the same parallelization
        fock_t* fockStateBuffer = 0;
        {
            fock_t* allowedFockStates = 0;
            //  Allocate memory to store part of the Fock basis on each node, and
            //  on the master node allocate enough space to store the full basis
            allowedFockStates = new (std::nothrow) fock_t[sectorCounter];
            fock_t* p_list = allowedFockStates;
            fock_t currState = utilities::binary::GenerateHammingNumber(nbrParticles, mpi.m_firstTask);
            for(fock_t i=0; i<fullNodeDim; ++i)
            {
                if(Constraint(currState))
                {
                    *p_list = currState;
                    ++p_list;
                }
                currState = utilities::binary::NextHammingNumber64(currState);
            }
            //  MPI gather the full list on the master node
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                fockStateBuffer = new fock_t[sectorTotalDimension];
            }
            MPI_Status status;
            mpi.Gather<fock_t>(allowedFockStates, sectorCounter, fockStateBuffer, sectorTotalDimension,
                               0, mpi.m_comm, status);
            delete[] allowedFockStates;
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
            m_basis.resize(cumulativeDim);
            //  Redistribute the fockStateBuffer over the other nodes
            fock_t taskList[mpi.m_nbrProcs];
            fock_t nodeDimList[mpi.m_nbrProcs];
            MPI_Gather(&m_firstState, 1, mpi.GetType<fock_t>(), taskList, 1, mpi.GetType<fock_t>(), 0, mpi.m_comm);
            MPI_Gather(&cumulativeDim, 1, mpi.GetType<fock_t>(), nodeDimList, 1, mpi.GetType<fock_t>(), 0, mpi.m_comm);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                for(int i=1; i<mpi.m_nbrProcs; ++i)
                {
                    MPI_Send(fockStateBuffer+taskList[i], nodeDimList[i], mpi.GetType<fock_t>(), i, 20+i, mpi.m_comm);
                }
                memcpy(m_basis.data(), fockStateBuffer, nodeDimList[0]*sizeof(fock_t));
                delete[] fockStateBuffer;
            }
            else
            {
                MPI_Recv(&m_basis[0], cumulativeDim, mpi.GetType<fock_t>(), 0, 20+mpi.m_id, mpi.m_comm, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(mpi.m_comm);
        m_basisStored = true;
	    return false;
    }
    
    //!
    //! Return a copy of the Fock basis to the given buffer
    //!
    void FermionFockBasis::GetFockBasis(
        fock_t* buffer,                 //!<    Buffer to copy to
        const fock_t dim,               //!<    Dimension of the given buffer
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        const
    {
        if(m_basis.size()>=dim)
        {
            memcpy(buffer, m_basis.data(), dim*sizeof(fock_t));
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
    void FermionFockBasis::SetFockBasis(
        fock_t* buffer,             //!<    Buffer to copy from
        const fock_t dim,           //!<    Dimension of the given buffer
        const fock_t firstState)    //!<    Index of first state represented
    {
        m_basis.resize(dim);
        memcpy(m_basis.data(), buffer, dim*sizeof(fock_t));
        //  Enforce ascending Fock state order
        std::sort(m_basis.begin(), m_basis.end());
        m_firstState = firstState;
        m_basisStored = true;
    }
}   //  End namespace diagonalization    
