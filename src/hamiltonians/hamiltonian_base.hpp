////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a a base class implementation to contain common 
//!     functions used in Hamiltonian type classes
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

#ifndef _HAMILTONIAN_BASE_HPP_INCLUDED_
#define _HAMILTONIAN_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "hamiltonian_data.hpp"
#include "../utilities/wrappers/mpi_wrapper.hpp"
#include "../utilities/data_structures/crs_sparse_matrix.hpp"
#include "../utilities/data_structures/mapped_sparse_matrix.hpp"                                               
#include "../utilities/mathematics/sparse_linear_algebra.hpp"
#include "../utilities/wrappers/arpack_wrapper.hpp" 
#include "../utilities/general/timer.hpp"
#include "../utilities/general/load_bar.hpp"
#include "../utilities/general/cout_tools.hpp"
#include "../utilities/wrappers/lapack_wrapper.hpp"
#include "../utilities/general/template_tools.hpp"
#include "../utilities/general/run_script.hpp"
#include "../utilities/wrappers/io_wrapper.hpp"
#include <functional>
#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

namespace diagonalization
{
    //  Pre-declare the MatrixVectorFunction class (see matrix_vector_routines.h 
    //  for implementation)
    template <typename T,class L1,class L2>
    class MatrixVectorFunction;

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief The Hamiltonian class contains common utility functions for
    //! matrix diagonalization
    //!
    //! The Hamiltonian is constructed from the outer product:
    //!
    //! |leftVector><rightVector|
    //!
    //! Where <rightVector| is given by applying the Hamiltonian to the full  
    //! Fock space and |leftVector> is the full set of Fock states
    //!
    //! The template parameter T denotes the data type of the matrix elements
    //! The template parameter B denotes the data type for the Fock basis
    //!
    ////////////////////////////////////////////////////////////////////////////////
    template <typename T,class B>
    class HamiltonianBase
    {
        typedef utilities::MappedSparseMatrix<T> mappedMatrix_t;
        typedef utilities::CrsSparseMatrix<T> crsSparseMatrix_t;
        public:
        HamiltonianData m_data;             //!<    Store Hamiltonian parameters
        B m_fockBasis;                      //!<    Store a lexicographically ordered 
                                            //!     Fock basis
        std::vector<double> m_eigenvalues;  //!<    Eigenvalues
        std::vector<T> m_eigenvectors;      //!<    Eigenvectors
        std::vector<T> m_fullMatrix;        //!<    To store the full matrix
                                            //!     representation of the Hamiltonian    
        mappedMatrix_t m_sparseMap;         //!<    Mapped Sparse matrix storage 
                                            //!     for fast population of our Hamiltonian matrix
        crsSparseMatrix_t m_crsSparseMatrix;//!<    Compressed Row storage format 
                                            //!     for fast matrix-vector operations on our matrix
        utilities::linearAlgebra::ArpackWrapper m_arpackWrapper;
                                            //!<    Wrapper class for ARPACK options
        protected:

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Allocate memory to store the full matrix representation of the 
        //! Hamiltonian
        //!
        //! If called with _PARALLEL_ flag then MUST be called on all nodes
        ////////////////////////////////////////////////////////////////////////////////
        void AllocateMatrixMemory(
            const utilities::storageMethod_t storageMethod, 
                                                //!<    Type of memory to allocate
            const fock_t nbrNonZeros,           //!<    Allocate a specific number of non-zeros
                                                //!     per matrix row for sparse storage
            const parallelFlag_t flag,          //!<    Flag whether to allocate
                                                //!     on all nodes, or just the master node
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            if(!m_data.m_matrixAllocated)
	        {
	            //  When parallelized, the Hamiltonian matrix can be divided up into strips
                //  with each strip being only addressed on a given node
	            if(_PARALLEL_ == flag || (_SERIAL_ == flag && 0 == mpi.m_id))
                {
	                if(utilities::_DENSE_ == storageMethod)
	                {
                        if(0 == mpi.m_id)	// FOR THE MASTER NODE
                        {
                            utilities::cout.DebuggingInfo()<<"\n\t- ALLOCATING "<<(m_data.m_fockSpaceDim*m_data.m_fockSpaceDim*(sizeof(T))/(1024.0*1024.0));
                            utilities::cout.DebuggingInfo()<<" MB TO STORE FULL HAMILTONIAN"<<std::endl;
                            if(m_data.m_fockSpaceDim*m_data.m_fockSpaceDim*(sizeof(T))>_MEMORY_QUERY_LIMIT_)
                            {
                                utilities::cout.DebuggingInfo()<<"\n\tWARNING: EXCEEDING "<<_MEMORY_QUERY_LIMIT_/(1024*1024*1024)<<" GB";
                                utilities::cout.DebuggingInfo()<<" MEMORY LIMIT! PRESS ANY KEY TO CONTINUE"<<std::endl;
                                getchar();
                            }
                        }
	                    m_fullMatrix.resize(m_data.m_nodeDim*m_data.m_fockSpaceDim);
	                }
	                else if(utilities::_SPARSE_CRS_ == storageMethod)
	                {
                        m_crsSparseMatrix.Initialize(m_data.m_nodeDim, m_data.m_fockSpaceDim, nbrNonZeros);
                    }
                    else if(utilities::_SPARSE_MAPPED_ == storageMethod)
	                {
                        m_sparseMap.Initialize(m_data.m_nodeDim, m_data.m_fockSpaceDim, nbrNonZeros);
                    }
                }
                m_data.m_matrixAllocated = true;
            }
            if(_PARALLEL_ == flag)
            {  
                MPI_Barrier(mpi.m_comm);
            }
	        return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Deallocate memory to store the full matrix representation 
        //! of the Hamiltonian
        ////////////////////////////////////////////////////////////////////////////////
        void DeallocateMatrixMemory(
            const utilities::storageMethod_t storageMethod)        //!<    Type of memory to de-allocate
        {
            if(m_data.m_matrixAllocated)
	        {
                utilities::cout.DebuggingInfo()<<"\n\t- DEALLOCATING EXISTING MATRIX STORAGE"<<std::endl;
	            if(utilities::_DENSE_==m_data.m_storageMethod)
	            {	
		            m_fullMatrix.clear();
                }
                else if(utilities::_SPARSE_MAPPED_ == storageMethod)
                {
                    m_sparseMap.Clear();
                }
                else if(utilities::_SPARSE_CRS_ == storageMethod)
                {
                    m_crsSparseMatrix.Clear();
                }
	
                m_data.m_matrixAllocated = false; 
            }
	        return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Converts the storage method to a new one (with sparse or
        //! dense)
        //!
        //!	If the matrix has already been allocated, then we transfer the 
        //!	current stored matrix to the new memory format
        //!
        //! This function MUST be called on all nodes if _PARALLEL_ is set
        ////////////////////////////////////////////////////////////////////////////////
        void ResetMatrixStorageMethod(
            const utilities::storageMethod_t newStorageMethod, 
                                                //!<    New storage method
            const parallelFlag_t flag,          //!<    Flag whether the reset is for
                                                //!     parallel re-allocation or
                                                //!     only re-allocates on the master node
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            utilities::Timer timer;
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
                {
                    if(utilities::_DENSE_ == newStorageMethod)
                    {
                        utilities::cout.AdditionalInfo()<<"\n\t- RESETTING STORAGE METHOD FROM SPARSE-MAPPED TO DENSE"<<std::endl;
                    }
                    else if(utilities::_SPARSE_CRS_ == newStorageMethod)
                    {
                        utilities::cout.AdditionalInfo()<<"\n\t- RESETTING STORAGE METHOD FROM SPARSE-MAPPED TO SPARSE-CRS"<<std::endl;
                    }
                }
                else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
                {
                    if(utilities::_DENSE_ == newStorageMethod)
                    {
                        utilities::cout.AdditionalInfo()<<"\n\t- RESETTING STORAGE METHOD FROM SPARSE-CRS TO DENSE"<<std::endl;
                    }
                    else if(utilities::_SPARSE_MAPPED_ == newStorageMethod)
                    {
                        utilities::cout.AdditionalInfo()<<"\n\t- RESETTING STORAGE METHOD FROM SPARSE-CRS TO SPARSE-MAPPED"<<std::endl;
                    }
                }
                else if(utilities::_DENSE_ == m_data.m_storageMethod)
                {
                    if(utilities::_SPARSE_CRS_ == newStorageMethod)
                    {
                        utilities::cout.AdditionalInfo()<<"\n\t- RESETTING STORAGE METHOD FROM DENSE TO SPARSE_CRS"<<std::endl;
                    }
                    else if(utilities::_SPARSE_MAPPED_ == newStorageMethod)
                    {
                        utilities::cout.AdditionalInfo()<<"\n\t- RESETTING STORAGE METHOD FROM DENSE TO SPARSE-MAPPED"<<std::endl;
                    }
                }
                timer.Start();
            }
            if(_PARALLEL_ == flag)
            {
                MPI_Barrier(mpi.m_comm);
            }
            if(m_data.m_matrixAllocated)
            {
	            if(newStorageMethod != m_data.m_storageMethod)	//	new method!
	            {
	                //  Enable memory allocation for the new storage method
		            m_data.m_matrixAllocated = false;
		            //	Copy matrix value across
		            if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
                    {               			
                        if(utilities::_DENSE_ == newStorageMethod)
                        {
                            this->AllocateMatrixMemory(newStorageMethod, m_sparseMap.GetNbrNonZeros(), flag, mpi);
                            //  Call the MappedSparseMatrix class's conversion function
                            m_sparseMap.ConvertToArray(m_fullMatrix.data());
                        }
                        else if(utilities::_SPARSE_CRS_ == newStorageMethod)
                        {
                            //  Call the CrsSparseMatrix class's initialization from m_sparseMap
                            m_crsSparseMatrix.Initialize(&m_sparseMap);
                            m_data.m_matrixAllocated = true;
                        }
                    }
                    else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
                    {
                        if(utilities::_DENSE_ == newStorageMethod)
                        {
                            this->AllocateMatrixMemory(newStorageMethod, m_crsSparseMatrix.GetNbrNonZeros(), flag, mpi);
                            //  Call the CrsSparseMatrix class's conversion function
                            m_crsSparseMatrix.ConvertToArray(m_fullMatrix.data(), m_data.m_fockSpaceDim);
                        }
                        else if(utilities::_SPARSE_MAPPED_ == newStorageMethod)
                        {
                            //  Call the MappedSparseMatrix class's initialization from a CRS matrix
                            m_sparseMap.Initialize(&m_crsSparseMatrix, m_data.m_nodeDim, m_data.m_fockSpaceDim);
                            m_data.m_matrixAllocated = true;
                        }
                    }
                    else if(utilities::_DENSE_ == m_data.m_storageMethod)
                    {
                        if(utilities::_SPARSE_CRS_ == newStorageMethod)
                        {
                            //  Call the CrsSparseMatrix class's initialization from an array
                            m_crsSparseMatrix.Initialize(m_fullMatrix.data(), m_data.m_nodeDim, m_data.m_fockSpaceDim, _SPARSE_ZERO_TOL_); 
                            m_data.m_matrixAllocated = true;
                        }
                        else if(utilities::_SPARSE_MAPPED_ == newStorageMethod)
                        {
                            //  Call the MappedSparseMatrix class's initialization from an array
                            m_sparseMap.Initialize(m_fullMatrix.data(), m_data.m_nodeDim, m_data.m_fockSpaceDim, _SPARSE_ZERO_TOL_);
                            m_data.m_matrixAllocated = true;
                        }
                    }
                    if(_PARALLEL_ == flag)
                    {
                        MPI_Barrier(mpi.m_comm);
                    }
		            //	Deallocate memory assigned to the old storage type
		            this->DeallocateMatrixMemory(m_data.m_storageMethod);
		            //  Update storage method flags
		            m_data.m_storageMethod = newStorageMethod;
		            m_data.m_matrixAllocated = true;
	            }
            }
            else
            {
                //	Update the storage format flag
                m_data.m_storageMethod = newStorageMethod;
            }
            if(0 == mpi.m_id)
            {
                timer.Stop();
            }
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Combine partial Hamiltonians stored on each separate node into 
        //!  a single combined Hamiltonian
        //!
        //! NOTE: There is no merge function for a hash map, so this type
        //! will first be converted to a CRS type
        //!
        //! This function MUST be called on all nodes if run in parallel
        ////////////////////////////////////////////////////////////////////////////////
        void MergeDistributedHamiltonian(
            utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            if(m_data.m_matrixAllocated && mpi.m_nbrProcs > 1)
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\n\t- GATHERING DISTRIBUTED HAMILTONIAN"<<std::endl;
                }
                if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
                {
                    //  Reset the storage method to CCS because this will allow simple MPI 
                    //  transfer of the underlying index and value arrays
                    this->ResetMatrixStorageMethod(utilities::_SPARSE_CRS_, _PARALLEL_, mpi);
                }
                //  Now gather the matrix using the SPARSE CRS or DENSE format
                if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
                {
                    //  Merge sparse matrix allocation from each node
                    m_crsSparseMatrix.MpiGather(0, mpi);
                }
                else if(utilities::_DENSE_ == m_data.m_storageMethod)
                {
                    //  Merge dense matrix allocation from each node
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {
                        utilities::cout.DebuggingInfo()<<"\n\t- ALLOCATING "<<((fock_t)m_data.m_fockSpaceDim*m_data.m_fockSpaceDim*(sizeof(T))/(1024.0*1024.0));
	                    utilities::cout.DebuggingInfo()<<" MB TO STORE BUFFER FOR MPI GATHER"<<std::endl;
	
	                    if((fock_t)(m_data.m_fockSpaceDim*m_data.m_fockSpaceDim*(sizeof(T)))>_MEMORY_QUERY_LIMIT_)
	                    {
	                        std::cerr<<"\n\tWARNING: EXCEEDING "<<_MEMORY_QUERY_LIMIT_/(1024*1024*1024)<<" GB";
	                        std::cerr<<" MEMORY LIMIT! PRESS ANY KEY TO CONTINUE"<<std::endl;
	                        getchar();
	                    }
                        //  Reallocate the m_fullMatrix address on the master node
                        //  to store the full matrix
                        m_fullMatrix.resize(m_data.m_fockSpaceDim*m_data.m_fockSpaceDim);
                    }
                    MPI_Barrier(mpi.m_comm);
                    MPI_Status status;
                    mpi.Gather<T>(m_fullMatrix.data(), m_data.m_nodeDim*m_data.m_fockSpaceDim, m_fullMatrix.data(), 
                                  m_data.m_fockSpaceDim*m_data.m_fockSpaceDim, 0, mpi.m_comm, status);
                    MPI_Barrier(mpi.m_comm);
                }
                if(0 != mpi.m_id)
                {
                    this->DeallocateMatrixMemory(m_data.m_storageMethod);
                }
            }
            //  Reset the node dimension data to describe a matrix stored on 
            //  a single node only
            m_data.SetDataDimensions(m_data.m_fockSpaceDim, 1, mpi);
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Run a check to see if the matrix is diagonal. LAPACK has trouble
        //! processing diagonal matrices, but then if the matrix is already
        //! diagonal we know what its eigenspectrum is anyway. 
        //!
        //! This function MUST be called on all nodes if run in parallel
        ////////////////////////////////////////////////////////////////////////////////
        bool MatrixDiagonalCheck(
            const parallelFlag_t flag,  //!<    Flag to set parallel implementation or not
            utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
        {
            int falseCount = 0;
            int totalFalseCount = 0;
            if(utilities::_DENSE_ == m_data.m_storageMethod && _SERIAL_ == flag)
            {
                for(fock_t i=0; i<m_data.m_fockSpaceDim; ++i)
                {
                    for(fock_t j=i+1; j<m_data.m_fockSpaceDim; ++j)
                    {
                        if(utilities::is_same<T,dcmplx>::value)
                        {
                            if(abs(m_fullMatrix[i*m_data.m_fockSpaceDim+j])>_SPARSE_ZERO_TOL_)
                            {
                                return false;
                            }
                        }
                        else if(utilities::is_same<T, double>::value)
                        {
                            if(fabs(m_fullMatrix[i*m_data.m_fockSpaceDim+j])>_SPARSE_ZERO_TOL_)
                            {
                                return false;
                            }
                        }
                    }
                }
                return true;
            }
            else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod && _PARALLEL_ == flag)
            {
                //  Determine the starting matrix index on each node
                mpi.DivideTasks(mpi.m_id, m_data.m_fockSpaceDim, mpi.m_nbrProcs,
                                &mpi.m_firstTask, &mpi.m_lastTask, false);
                //  Iterate over all off diagonal sparse matrix elements  
                int* p_rowStarts = m_crsSparseMatrix.GetRowStarts();
                int* p_cols      = m_crsSparseMatrix.GetCols();
                T* p_data        = m_crsSparseMatrix.GetData();
                for(int i=0; i<m_crsSparseMatrix.GetNbrRows(); ++i)
                {
                    for(int j = p_rowStarts[i]; j<p_rowStarts[i+1]; ++j)
                    {
                        if(utilities::is_same<T,dcmplx>::value)
                        {
                            if((p_cols[j]!=(int)mpi.m_firstTask+i) && abs(p_data[j])>_SPARSE_ZERO_TOL_)
                            {
                                falseCount = 1;
                                goto end;
                            }
                       }
                       else if(utilities::is_same<T,double>::value)
                       {
                            if((p_cols[j]!=(int)mpi.m_firstTask+i) && fabs(p_data[j])>_SPARSE_ZERO_TOL_)
                            {
                                falseCount = 1;
                                goto end;
                            }
                       }
                    }
                }
            }
            else
            {
                return false;
            }
            end:
            //  Count the number of false values set on all nodes
            MPI_Allreduce(&falseCount, &totalFalseCount, 1, mpi.GetType<int>(), MPI_SUM, mpi.m_comm);
            //  The matrix is diagonal if there were no false returns from any node
            return totalFalseCount == 0;
        }

        //////      Public interface functions      ////////////////////////////////////////////////////
	    public:

        //!
        //! Default constructor declaration
        //!
        HamiltonianBase()
        {}

        //!
        //! Copy constructor declaration
        //!
        HamiltonianBase(const HamiltonianBase& other)
        :
        m_data(other.m_data),
        m_fockBasis(other.m_fockBasis),
        m_fullMatrix(other.m_fullMatrix),
        m_sparseMap(other.m_sparseMap),
        m_crsSparseMatrix(other.m_crsSparseMatrix),
        m_eigenvalues(other.m_eigenvalues),
        m_eigenvectors(other.m_eigenvectors),
        m_arpackWrapper(other.m_arpackWrapper)
        {}

        //!
        //! Copy assignment operator declaration
        //!
        HamiltonianBase& operator=(const HamiltonianBase& other)
        {
            m_data          = other.m_data;
            m_fockBasis     = other.m_fockBasis;
            m_fullMatrix    = other.m_fullMatrix;
            m_sparseMap     = other.m_sparseMap;
            m_crsSparseMatrix = other.m_crsSparseMatrix;
            m_eigenvalues   = other.m_eigenvalues;
            m_eigenvectors  = other.m_eigenvectors;
            m_arpackWrapper = other.m_arpackWrapper;
            return *this;
        }

        //!
        //! Destructor declaration
        //! 
        ~HamiltonianBase()
        {}

        protected:

        //!
        //! Initialization function to override default construction
        //!
        void BaseInitialize(
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
                {
                    utilities::cout.SecondaryOutput()<<"\t\tMATRIX STORAGE METHOD:\tSPARSE MAPPED MATRIX"<<std::endl;
                }
                else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
                {
                    utilities::cout.SecondaryOutput()<<"\t\tMATRIX STORAGE METHOD:\tSPARSE COMPRESSED ROW STORAGE"<<std::endl;
                }
                else if(utilities::_DENSE_ == m_data.m_storageMethod)
                {
                    utilities::cout.SecondaryOutput()<<"\t\tMATRIX STORAGE METHOD:\tDENSE"<<std::endl;
                }
	            #if _USING_ITERATIVE_HAMMING_WEIGHT_ALGORITHM_
	                utilities::cout.AdditionalInfo()<<"\n\t\tUSING ITERATIVE HAMMING WEIGHT ALGORITHM"<<std::endl;
	            #else
	                utilities::cout.AdditionalInfo()<<"\n\t\tUSING ARITHMETIC HAMMING WEIGHT ALGORITHM"<<std::endl;
	            #endif
	            if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
	            {
	                #if _SPEED_OPTIMIZED_MAP_
	                    utilities::cout.AdditionalInfo()<<"\n\t\tUSING SPEED OPTIMIZED HASH MAP (google::dense_hash_map)"<<std::endl;
	                #elif _MEMORY_OPTIMIZED_MAP_
	                    utilities::cout.AdditionalInfo()<<"\n\t\tUSING MEMORY OPTIMIZED HASH MAP (google::sparse_hash_map)"<<std::endl; 
	               #endif        
	            }
	        }
        }

        //!
        //! Clear all class data structures
        //! 
        void BaseClear()
        {
            m_fockBasis.Clear();
            m_data.Clear();
            m_fullMatrix.clear();
            m_crsSparseMatrix.Clear();       
            m_sparseMap.Clear();
            m_eigenvalues.clear();
            m_eigenvectors.clear();
        }
        
        public:

        //!
        //! Initialize ARPACK options
        //! 
        void InitializeArpack(
            boost::program_options::variables_map* optionList,    //!<    Parsed command line argument list
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {
            if(0 == mpi.m_id)
            {
                m_arpackWrapper.SetEigenvalueMagnitude('S');
                if((*optionList)["use-initial"].as<bool>())
                {
                    m_arpackWrapper.UseInitialVector();
                }
                if((*optionList)["store-final"].as<bool>())
                {
                    m_arpackWrapper.StoreFinalVector();
                }
            }
            m_arpackWrapper.MpiSync(0, mpi);
        }
        
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Update ARPACK input and output file file names in order to change 
        //! the labels at any point
        ////////////////////////////////////////////////////////////////////////////////
        void UpdateArpackFileNames(
            const std::string inFileName,     //!<    Name of ARPACK vector input file
            const std::string outFileName,    //!<    Name of ARPACK vector output file
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        {          
            m_arpackWrapper.SetInitialVectorFile(inFileName);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<"\n\t SET ARPACK INITIAL VECTOR FILE: "<<inFileName<<std::endl;
            }
            m_arpackWrapper.SetFinalVectorFile(outFileName); 
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<"\n\t SET ARPACK FINAL VECTOR FILE: "<<outFileName<<std::endl;
            }
        }
  
        //!
        //! Set the Fock basis
        //!
        void SetFockBasis(const B& fockBasis)
        {
            m_fockBasis = fockBasis;
        }

        //!
        //! Set the number of eigenvalues parameter
        //!
        void SetNbrEigenvalues(
            const iSize_t nbrEigenvalues,           //!<    NUmber of eigenvalues to set
            const utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
        {
            if(nbrEigenvalues>m_data.m_fockSpaceDim)
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::cerr<<"\n\tWARNING in SpinlessFermionHamiltonian: nbrEigenvalues larger than Fock space dimension:"<<std::endl;
                    std::cerr<<"\tSet number of eigenvalues to "<<m_data.m_fockSpaceDim<<std::endl;
                }
                m_data.m_nbrEigenvalues = m_data.m_fockSpaceDim;
            }
            else
            {
                m_data.m_nbrEigenvalues = nbrEigenvalues;
            }
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Call a dense matrix diagonalization routine from LAPACK
        //!
        //! This function MUST be called on all nodes if run in parallel
        //!
        //! NOTE: deallocates matrix storage memory once completed
        //!
        //! \return true if success, false if failure
        ////////////////////////////////////////////////////////////////////////////////
        bool FullDiagonalize(
            iSize_t& nbrEigenvalues,        //!<    Number of eigenvalues to find
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        {
            MPI_Barrier(mpi.m_comm);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                 utilities::cout.SecondaryOutput()<<"\n\t============ DIAGONALIZING HAMILTONIAN WITH LAPACK ============ "<<std::endl;
            }
            //  Make sure that all nodes know how many eigenvalues we are including
            this->SetNbrEigenvalues(nbrEigenvalues, mpi);
            //  LAPACK will only work on a single node
            //  (use SCALAPACK for parallel algorithms?)
            this->MergeDistributedHamiltonian(mpi);
            if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod || utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
            {
                //  Need to convert to dense storage format first
                this->ResetMatrixStorageMethod(utilities::_DENSE_, _SERIAL_, mpi);
            }
            //  Proceed on the master node
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                if(m_data.m_matrixAllocated)
                {
                    m_eigenvalues.resize(m_data.m_nbrEigenvalues);
                    m_eigenvectors.resize(m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim);
                    if(this->MatrixDiagonalCheck(_SERIAL_,mpi))
                    {
                        utilities::cout.SecondaryOutput()<<"\n\tMATRIX WAS ALREADY DIAGONAL - GENERATING SPECTRUM DIRECTLY"<<std::endl;
                        //  For a diagonal matrix, the eigenspectrum can be easily extracted
                        double* diagBuffer = new double[m_data.m_fockSpaceDim];
                        double* sortBuffer = new double[m_data.m_fockSpaceDim];
                        //  Populate a buffer with diagonal entries
                        for(fock_t i=0; i<m_data.m_fockSpaceDim; ++i)
                        {
                            diagBuffer[i] = std::real(m_fullMatrix[i*m_data.m_fockSpaceDim+i]);
                        }
                        memcpy(sortBuffer, diagBuffer, m_data.m_fockSpaceDim*sizeof(double));
                        //  Put the values in descending order
                        std::sort(sortBuffer,sortBuffer+m_data.m_fockSpaceDim);
                        //  Read out the lowest eigenvalues
                        //  and generate their corresponding eigenvectors
                        //  by generating zero vectors and specifying
                        //  a 1.0 value for each degenerate value,
                        //  the normalizing at the end
                        for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                        {
                            m_eigenvalues[i] = sortBuffer[i];
                            double norm = 0.0;
                            for(fock_t j=0; j<m_data.m_fockSpaceDim; ++j)
                            {    
                                if(fabs(m_eigenvalues[i]-diagBuffer[j])<_SPARSE_ZERO_TOL_)
                                {
                                    m_eigenvectors[i*m_data.m_fockSpaceDim+j] = 1.0;
                                    norm += 1.0;
                                }
                                else
                                {
                                    m_eigenvectors[i*m_data.m_fockSpaceDim+j] = 0.0;
                                }
                            }
                            for(fock_t j=0; j<m_data.m_fockSpaceDim; ++j)
                            {
                                m_eigenvectors[i*m_data.m_fockSpaceDim+j] /= norm;
                            }
                        }
                        delete[] diagBuffer;
                        delete[] sortBuffer;
                    }
                    else
                    {
                        utilities::cout.SecondaryOutput()<<"\n\tCALLING LAPACK ROUTINE"<<std::endl;
                        utilities::Timer timer;     //  Time the operation
                        timer.Start();
                        utilities::linearAlgebra::DiagonalizeSymmetricMatrix<T>(m_fullMatrix.data(), m_eigenvectors.data(), 
                                                                                m_eigenvalues.data(), m_data.m_fockSpaceDim, 
                                                                                m_data.m_nbrEigenvalues, 'U');
                        timer.Stop();
	                }
	                //  Now deallocate the memory assigned to store the Hamiltonian
	                this->DeallocateMatrixMemory(m_data.m_storageMethod);
                    utilities::cout.AdditionalInfo()<<"\n\t============ HAMILTONIAN SUCESSFULLY DIAGONALIZED ============ "<<std::endl;
                    m_data.m_eigenvaluesCalculated  = true;
                    m_data.m_eigenvectorsCalculated = true;
                    this->PrintEigensystem();
                }
                else
                {
                    m_data.m_eigenvaluesCalculated  = false;
                    m_data.m_eigenvectorsCalculated = false;
                }
            }
            mpi.Sync(&m_data.m_eigenvaluesCalculated, 1, 0);
            mpi.Sync(&m_data.m_eigenvectorsCalculated, 1, 0);
            //  We need to scatter the eigenvectors over multiple nodes if we want to perform
            //  any further data analysis in parallel
            if(m_data.m_eigenvaluesCalculated)
            {
                //  Synchronize with the master node
                mpi.Sync(&m_eigenvalues, 0);
            }
            if(m_data.m_eigenvectorsCalculated)
            {
                //  Call MPI scatter on the eigenvectors
                //  First copy the existing eigenvectors into the scatter buffer
                T* eigenvectorsBuffer = 0;
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    eigenvectorsBuffer = new T[m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim];
                    memcpy(eigenvectorsBuffer, m_eigenvectors.data(), m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim*sizeof(T));
                }   
                mpi.DivideTasks(mpi.m_id, m_data.m_fockSpaceDim, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
                fock_t eigenvectorNodeDim = mpi.m_lastTask - mpi.m_firstTask + 1;
                m_eigenvectors.resize(m_data.m_nbrEigenvalues*eigenvectorNodeDim);
                //  MPI scatter each eigenvector
                for(iSize_t e=0; e<m_data.m_nbrEigenvalues; ++e)
                {
                    MPI_Status status;
                    mpi.Scatter<T>(eigenvectorsBuffer+e*m_data.m_fockSpaceDim, m_data.m_fockSpaceDim,
                                   m_eigenvectors.data()+e*eigenvectorNodeDim, eigenvectorNodeDim, 0, mpi.m_comm, status);
                }
                if(0!=eigenvectorsBuffer)   delete[] eigenvectorsBuffer;
            }
            return m_data.m_eigenvaluesCalculated && m_data.m_eigenvectorsCalculated;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Call a sparse Lanczos algorithm from the ARPACK/PARPACK library
        //!
        //! This function MUST be called on all nodes if run in parallel
        //!
        //! NOTE: deallocates matrix storage memory once completed
        //!
        //! \return true if success, false if failure
        ////////////////////////////////////////////////////////////////////////////////
        template <class L1, class L2>
        bool LanczosDiagonalize(
            L1* quadraticLookUpTables,          //!<    Class container for look-up tables
            L2* quarticLookUpTables,            //!<    Class container for look-up tables
            MatrixVectorFunction<T, L1, L2>& mv,//!<    Class to perform matrix-vector operation
            iSize_t& nbrEigenvalues,            //!<    Number of eigenvalues to find
            utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
        {
            ////////////////////////////////////////////////////////////////////////////////
            #if _ENABLE_ARPACK_
            MPI_Barrier(mpi.m_comm);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                 utilities::cout.SecondaryOutput()<<"\n\t============ DIAGONALIZING HAMILTONIAN WITH ARPACK ============ "<<std::endl;
            }
            //  Make sure that all nodes know how many eigenvalues we are including
            this->SetNbrEigenvalues(nbrEigenvalues,mpi);
            if(utilities::_DENSE_ == m_data.m_storageMethod || utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
            {
                //  Need to convert to sparse-ccs storage format first
                //  for fast matrix-vector operations and for the
                //  MergeDistributedHamiltonian function
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    utilities::cout.SecondaryOutput()<<"\n\t- RESETTING STORAGE METHOD TO CRS"<<std::endl;
                }
                this->ResetMatrixStorageMethod(utilities::_SPARSE_CRS_, _PARALLEL_, mpi);
            }
            if(m_data.m_matrixAllocated)
            {
                //  Allocate memory to store eigenvalues/eigenvectors
                //  Distributed over a number of different nodes
                m_eigenvalues.resize(m_data.m_nbrEigenvalues);
                m_eigenvectors.resize(m_data.m_nbrEigenvalues*m_data.m_nodeDim);
                if(this->MatrixDiagonalCheck(_PARALLEL_, mpi))
                {
                    if(0 == mpi.m_id)    //  For the master node
                    {
                        utilities::cout.SecondaryOutput()<<"\n\tMATRIX WAS ALREADY DIAGONAL - GENERATING SPECTRUM DIRECTLY"<<std::endl;
                    }
                    //  For a diagonal matrix, the eigenspectrum can be easily extracted
                    double* diagBuffer = new double[m_data.m_nodeDim];
                    //  Populate a buffer with diagonal entries
                    mpi.DivideTasks(mpi.m_id, m_data.m_fockSpaceDim, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
                    for(fock_t i=mpi.m_firstTask; i<=mpi.m_lastTask; ++i)
                    {
                        diagBuffer[i-mpi.m_firstTask] = std::real(m_crsSparseMatrix.Value(i-mpi.m_firstTask, i));
                    }
                    //  Merge the diagonal elements on the master node
                    double* mergeBuffer=0;
                    if(0 == mpi.m_id)    //  For the master node
                    {
                        mergeBuffer = new double[m_data.m_fockSpaceDim];
                    }
                    MPI_Status status;
                    mpi.Gather<double>(diagBuffer, m_data.m_nodeDim, mergeBuffer, m_data.m_fockSpaceDim, 0, mpi.m_comm, status);
                    if(0 == mpi.m_id)    //  proceed on the master node only
                    {
                        double* sortBuffer = new double[m_data.m_fockSpaceDim];
                        memcpy(sortBuffer, mergeBuffer, m_data.m_fockSpaceDim*sizeof(double));
                        //  Put the values in descending order
                        std::sort(sortBuffer, sortBuffer+m_data.m_fockSpaceDim);
                        //  Read out the lowest eigenvalues
                        //  and generate their corresponding eigenvectors
                        //  by generating zero vectors and specifying
                        //  a 1.0 value for each degenerate value,
                        //  the normalizing at the end
                        for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                        {
                            m_eigenvalues[i] = sortBuffer[i];
                        }
                        delete[] sortBuffer;
                        delete[] mergeBuffer;
                    }
                    //  Set the eigenvalues on all nodes
                    mpi.Sync(&m_eigenvalues, 0);
                    //  Set the eigenvectors
                    for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                    {        
                        double norm = 0.0;
                        double totalNorm = 0.0;
                        for(fock_t j=mpi.m_firstTask; j<=mpi.m_lastTask; ++j)
                        {    
                            if(fabs(m_eigenvalues[i]-diagBuffer[j-mpi.m_firstTask])<_SPARSE_ZERO_TOL_)
                            {
                                m_eigenvectors[i*m_data.m_nodeDim+j-mpi.m_firstTask] = 1.0;
                                norm += 1.0;
                            }
                            else
                            {
                                m_eigenvectors[i*m_data.m_nodeDim+j-mpi.m_firstTask] = 0.0;
                            }
                        }
                        //  Calculate the total norm
                        MPI_Allreduce(&norm, &totalNorm, 1, mpi.GetType<double>(), MPI_SUM, mpi.m_comm);
                        for(fock_t j=mpi.m_firstTask; j<=mpi.m_lastTask; ++j)
                        {
                            m_eigenvectors[i*m_data.m_nodeDim+j-mpi.m_firstTask] /= totalNorm;
                        }
                    }
                    MPI_Barrier(mpi.m_comm);
                    //  Clear buffer allocation
                    delete[] diagBuffer;
                }
                else
                {
                    //////      DIAGONALIZE WITH A LANCZOS ALGORITHM        ////////
                    //  Run some checks for known badly conditioned matrices to pass to ARPACK
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {   
                        utilities::cout.SecondaryOutput()<<"\n\t- CALLING (P)ARPACK ROUTINE"<<std::endl;
                        if(m_data.m_nodeDim<=3)
                        {
                            std::cerr<<"\n\tWARNING - CALLING (P)ARPACK FOR A VERY SMALL MATRIX MAY LEAD TO AN ERROR. USE LAPACK INSTEAD?"<<std::endl;
                        }
                    }
                    utilities::Timer timer;
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {
                        timer.Start();  //  Time the operation
                    }
                    //////  Prepare the CRS matrix format for parallel matrix-vector
                    //////  multiplication routine
                    m_crsSparseMatrix.PrepareForParallelMultiplication(&(mv.dataDistribution), &(mv.commGroups), 
                                                                       &(mv.inputVectorBuffer), &(mv.outputVectorBuffer), mpi);
                    //////  CALL MAIN ARPACK ROUTINE
                    m_arpackWrapper.DiagonalizeSymmetricMatrix<T>(std::bind(&MatrixVectorFunction<T, L1, L2>::Simple, 
                                                                  mv, std::placeholders::_1, std::placeholders::_2, 
                                                                  std::placeholders::_3), m_eigenvectors.data(), 
                                                                  m_eigenvalues.data(), m_data.m_nodeDim, m_data.m_fockSpaceDim,
                                                                  m_data.m_nbrEigenvalues, mpi);
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {
                        timer.Stop();
                    }
                }
                this->DeallocateMatrixMemory(m_data.m_storageMethod);
                MPI_Barrier(mpi.m_comm);
                if(0 == mpi.m_id)
                {   
                    utilities::cout.AdditionalInfo()<<"\n\t============ HAMILTONIAN SUCESSFULLY DIAGONALIZED ============ "<<std::endl;
                }
                m_data.m_eigenvaluesCalculated  = true;
                m_data.m_eigenvectorsCalculated = true;
                if(0 == mpi.m_id)
                {
                    this->PrintEigensystem();
                }
            }
            else
            {
	            //	ERROR
	            m_data.m_eigenvaluesCalculated  = false;
                m_data.m_eigenvectorsCalculated = false;
            }
            return m_data.m_eigenvaluesCalculated && m_data.m_eigenvectorsCalculated;
            #else
            //  If ARPACK is not enabled, default to return a flag to indicate that 
            //  diagonalization failed
            return false;
	        #endif
	        ////////////////////////////////////////////////////////////////////////////////
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief PrintHamiltonian out the matrix to screen
        //!
        //! Only PrintHamiltonian fully for small matrices
        ////////////////////////////////////////////////////////////////////////////////
        void PrintHamiltonian(
            const int nodeId,                   //!<    Node to print on
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
            const
        {
            const int fullPrintLimit     = 20;
            const int nonZerosPrintLimit = 60;
        
            if(nodeId == mpi.m_id && m_data.m_matrixAllocated)
            {
	            if(m_data.m_fockSpaceDim<fullPrintLimit)   
	            {
		            utilities::cout.SecondaryOutput()<<"\n\t- PRINT HAMILTONIAN"<<std::endl<<std::endl;
		            for(fock_t i=0; i<m_data.m_nodeDim; ++i)
		            {
			            utilities::cout.SecondaryOutput()<<"\tROW "<<i+mpi.m_firstTask<<"\t";
			            for(fock_t j=0; j<m_data.m_fockSpaceDim; ++j)
			            {
			                if(utilities::_DENSE_ == m_data.m_storageMethod)
                            {   
			                    utilities::cout.SecondaryOutput()<<m_fullMatrix[i*m_data.m_fockSpaceDim+j]<<" ";
			                }
			                else if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
			                {
			                    utilities::cout.SecondaryOutput()<<m_sparseMap.Value(i,j)<<" ";
			                }
			                else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
			                {
			                    utilities::cout.SecondaryOutput()<<m_crsSparseMatrix.Value(i,j)<<" ";
			                }
			            }
			            utilities::cout.SecondaryOutput()<<std::endl;
		            }
		        }
		        else if(m_data.m_fockSpaceDim<nonZerosPrintLimit)  //  PrintHamiltonian only non-zero elements
                {
                    utilities::cout.SecondaryOutput()<<"\n\t- PRINT HAMILTONIAN NON-ZEROS (ROW,COL) VALUE"<<std::endl<<std::endl;
                    if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
			        {
			            m_sparseMap.Print();
                    }
                    else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
			        {
			            m_crsSparseMatrix.Print();
			        }
                }
            }
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief PrintHamiltonian the amount of memory allocated to store the matrix
        //!
        //! This function MUST be called on all nodes if run in parallel
        ////////////////////////////////////////////////////////////////////////////////
        void PrintMemoryAllocation(
            const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class#
            const
        {
            MPI_Barrier(mpi.m_comm);
            if(m_data.m_matrixAllocated)
            {
                if(utilities::_SPARSE_MAPPED_ == m_data.m_storageMethod)
                {
                    utilities::cout.DebuggingInfo()<<"\n\t- ON NODE "<<mpi.m_id<<" CURRENTLY HAVE ALLOCATED "<<(m_sparseMap.GetNbrNonZeros()*(m_sparseMap.SizeOfElement())/(1024.0*1024.0));
                    utilities::cout.DebuggingInfo()<<" MB TO STORE SPARSE HAMILTONIAN ("<<m_sparseMap.GetNbrNonZeros()<<" NON-ZERO ELEMENTS)."<<std::endl;
                }
                else if(utilities::_SPARSE_CRS_ == m_data.m_storageMethod)
                {
                    utilities::cout.DebuggingInfo()<<"\n\t- ON NODE "<<mpi.m_id<<" CURRENTLY HAVE ALLOCATED "<<(m_crsSparseMatrix.GetNbrNonZeros()*(sizeof(T)+sizeof(crsIndex_t))/(1024.0*1024.0));
                    utilities::cout.DebuggingInfo()<<" MB TO STORE SPARSE HAMILTONIAN ("<<m_crsSparseMatrix.GetNbrNonZeros()<<" NON-ZERO ELEMENTS)."<<std::endl;
                }
                else if(utilities::_DENSE_ == m_data.m_storageMethod)
                {
                    utilities::cout.DebuggingInfo()<<"\n\t- ON NODE "<<mpi.m_id<<" CURRENTLY HAVE ALLOCATED "<<(m_data.m_fockSpaceDim*m_data.m_fockSpaceDim*(sizeof(T))/(1024.0*1024.0));
                    utilities::cout.DebuggingInfo()<<" MB TO STORE DENSE HAMILTONIAN"<<std::endl;
                }
            }
            MPI_Barrier(mpi.m_comm);
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to print the eigenvalues and eigenvectors
        //!
        ////////////////////////////////////////////////////////////////////////////////
        void PrintEigensystem() const
        {
            if(m_data.m_eigenvectorsCalculated)
            {
                utilities::cout.MainOutput()<<"\n\tEIGENVALUES:\n"<<std::endl;
                for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                {   
                    utilities::cout.MainOutput().precision(10);
                    utilities::cout.MainOutput()<<"\t\t"<<m_eigenvalues[i]<<std::endl;
                }
            }
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Return the currently stored node dimension of the eigenvalues buffer
        //!
        ////////////////////////////////////////////////////////////////////////////////
        fock_t GetEigenvectorNodeDim() const
        {
            return m_eigenvectors.size()/m_data.m_nbrEigenvalues;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Generate a plot of the matrix
        //!
        //! This function can be called in parallel it only makes a plot of the
        //! part of the Hamiltonian currently stored on the nodes it runs on
        //!
        //! \return true if the calculation was successful, false otherwise
        ////////////////////////////////////////////////////////////////////////////////
        bool PlotHamiltonian(
            const std::string fileName,     //!<    Name of the file where Hamiltonian plot is stored
            const std::string tmpDir,       //!<    A temporary directory where data are stored
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        {
            if(m_data.m_matrixAllocated)
            {
                if(0 == mpi.m_id) //  FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\n\t- PLOT HAMILTONIAN"<<std::endl;
                }
	            //  Output the Hamiltonian into a file
	            std::stringstream fileNameStream;
	            fileNameStream.str("");
	            fileNameStream<<tmpDir<<"hamiltonianPlotData_id_"<<mpi.m_id<<".tmp";
	            io::fileFormat_t format = io::_TEXT_;
                this->HamiltonianToFile(fileNameStream.str(), format, mpi);
	            //  Generate the python script that will perform the plotting
	            std::stringstream pythonScript;
	            pythonScript.str();
	            pythonScript<<"#! //usr//bin//env python"<<_PYTHON_VERSION_<<"\n"
	            "import matplotlib                              \n"
	            "import numpy as np                             \n"
	            "import math                                    \n"
	            "matplotlib.use('Agg')                        \n\n"
	            "import matplotlib.pyplot as plt              \n\n"
	            "data = []                                      \n"
                "fin=open(\""<<fileNameStream.str()<<"\")     \n\n"
	            "dimX = float(fin.readline())                   \n"
	            "dimY = float(fin.readline())                 \n\n"
	            "for line in fin:                               \n"
	            "    data.append(float(line))                 \n\n"
	            "fin.close()                                  \n\n"
	            "data2 = np.asarray(data).reshape(dimX,dimY)  \n\n"
	            "plt.imshow(data2,interpolation='none')         \n"
	            "plt.colorbar()                                 \n"
	            "plt.title(\"Matrix plot of Hamiltonian\")\n"
	            "plt.savefig(\""<<fileName<<"id_"<<mpi.m_id<<".pdf\",bbox_inches=\'tight\') \n"
	            "plt.close() \n";
	            //  Execute the script
                utilities::Script myScript;
                myScript.SetScript(pythonScript.str());
                myScript.Execute();
	            if(0 == mpi.m_id) //  FOR THE MASTER NODE
                {
	                utilities::cout.AdditionalInfo()<<"\n\tHAMITLONIAN MATRIX PLOT(S) SUCESSFULLY GENERATED."<<std::endl;
	            }
            }
            else
            {
                if(0 == mpi.m_id) //  FOR THE MASTER NODE
                {
                    std::cerr<<"\n\t- PLOT HAMILTONIAN ERROR: not allocated/no longer allocted!"<<std::endl;
                }
                return false;
            }
            return true;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Put the Hamiltonian data in a file in a co-ordinate matrix format
        //!
        //! This function MUST be called on all nodes if run in parallel 
        //!
        //! NOTE: multiple files will contain the rows of the matrix
        //! stored on each node (make sure that different nodes use
        //! different file names)
        ////////////////////////////////////////////////////////////////////////////////
        void HamiltonianToFile(
            const std::string fileName,     //!<    Name of output file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        {
            if(m_data.m_matrixAllocated)
            {
                if(0 == mpi.m_id) //  FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\t- WRITING FILE"<<std::endl;
                }
                //  Need to convert to sparse CRS storage format first
                if(utilities::_SPARSE_CRS_ != m_data.m_storageMethod)
                {
                    this->ResetMatrixStorageMethod(utilities::_SPARSE_CRS_, _PARALLEL_, mpi);
                }
                m_crsSparseMatrix.ToFile(fileName, format, mpi);
                if(0 == mpi.m_id) //  FOR THE MASTER NODE
                {
                    utilities::cout.AdditionalInfo()<<"\n\t- DONE"<<std::endl;
                }
            }
            else
            {
                if(0 == mpi.m_id) //  FOR THE MASTER NODE
                {
                    std::cerr<<"\n\t- Hamiltonain not allocated - failed to srite to file"<<std::endl;
                }
            }
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Read the Hamiltonian data from a file in a co-ordinate matrix format
        //!
        //! This function MUST be called on all nodes if run in parallel 
        //!
        //! NOTE: multiple files will contain the rows of the matrix
        //! stored on each node (make sure that different nodes use
        //! different file names)
        ////////////////////////////////////////////////////////////////////////////////
        void HamiltonianFromFile(
            const std::string fileName,     //!<    Name of input file
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        {
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<"\n\t- RETRIEVE HAMILTONIAN FROM DISK"<<std::endl;
            }
            //  File format is always CRS sparse storage
            m_data.m_storageMethod = utilities::_SPARSE_CRS_;
            m_crsSparseMatrix.FromFile(fileName, format, mpi);
            m_data.m_nodeDim = m_crsSparseMatrix.GetNbrRows();
            m_data.m_matrixAllocated = true;
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<"\n\t- DONE"<<std::endl;
            }
            return;
        }
        
        //!
        //! Get the address of the eigenvalue array
        //!
        void GetEigenvalues(
            double* buffer,         //!<    Buffer to store returned eigenvalues
            const iSize_t nbrEigenvalues)
                                    //!<    Number of eigenvalues
            const
        {
            if(nbrEigenvalues > m_eigenvalues.size())
            {
                std::cerr<<"\n\t- ERROR in GetEigenvalues: incorrect number specified"<<std::endl;
            }  
            memcpy(buffer,m_eigenvalues.data(), sizeof(double)*nbrEigenvalues);
        }

        //!
        //! Output all stored eigenvalue and eigenvector data to a text file
        //!
        void EigensystemToFile(
            const std::string fileName,     //!<    Name of output file
            bool writeValues,               //!<    Optionally write eigenvalues
            bool writeVectors,              //!<    Optionally write eigenvectors 
                                            //!     (can be turned off to save disk space)
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        {   
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<"\n\t- WRITING EIGENSYSTEM DATA TO FILE "<<fileName.c_str()<<std::endl;
            }
            //  Set flag to include eigenvalues with eigenvectors
            if(writeVectors)     
            {
                writeValues = writeVectors;
            }
            if(!m_data.m_eigenvaluesCalculated && writeValues)
            {
                std::cerr << "ERROR eigenvalues cannot be written to file as they are not calculated" << std::endl;
                mpi.m_exitFlag=true;
                return;
            }
            if(!m_data.m_eigenvectorsCalculated && writeVectors)
            {
                std::cerr << "ERROR eigenvectors cannot be written to file as they are not calculated" << std::endl;
                mpi.m_exitFlag=true;
                return;
            }
            //  Call MPI gather on the eigenvectors, if writing to file
            T* eigenvectorsBuffer = 0;
            if(writeVectors)
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    eigenvectorsBuffer = new T[m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim];
                }
                for(iSize_t e=0; e<m_data.m_nbrEigenvalues; ++e)
                {
                    MPI_Status status;
                    fock_t eigenvectorNodeDim = this->GetEigenvectorNodeDim();
                    mpi.Gather<T>(m_eigenvectors.data()+e*eigenvectorNodeDim, eigenvectorNodeDim,
                                  eigenvectorsBuffer+e*m_data.m_fockSpaceDim, m_data.m_fockSpaceDim, 0, mpi.m_comm, status);
                }
            }
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                std::ofstream f_out;
                utilities::GenFileStream(f_out, fileName, format, mpi);
                if(!mpi.m_exitFlag)
                {
                    if(writeValues)
                    {
                        if(io::_TEXT_ == format)
                        {
                            f_out << "##  NBR EIGENVALUES\n";
                            f_out << m_data.m_nbrEigenvalues << "\n";
                            f_out << "##  EIGENVALUES\n";
                            for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                            {
                                f_out << std::setprecision(15) << m_eigenvalues[i] << "\n";
                            }
                        }
                        else if(io::_BINARY_ == format)
                        {
                            f_out.write((char*)&(m_data.m_nbrEigenvalues), sizeof(unsigned int));
                            f_out.write((char*)m_eigenvalues.data(), m_data.m_nbrEigenvalues*sizeof(double));
                        }
                    }
                    if(writeVectors)
                    {
                        if(io::_TEXT_ == format)
                        {
                            f_out << "##  FOCK SPACE DIMENSION\n";
                            f_out << m_data.m_fockSpaceDim << "\n";
                            f_out << "##  EIGENVECTORS\n";
                            for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                            {
                                T* p_eigs = eigenvectorsBuffer + i*m_data.m_fockSpaceDim;
                                for(fock_t j=0; j<m_data.m_fockSpaceDim; ++j, ++p_eigs)
                                {
                                    std::string streamData = utilities::ToStream(*p_eigs);
                                    f_out << streamData << "\n";
                                }
                            }
                        }
                        else if(io::_BINARY_ == format)
                        {
                            f_out.write((char*)&(m_data.m_fockSpaceDim), sizeof(unsigned int));
                            f_out.write((char*)eigenvectorsBuffer, m_data.m_fockSpaceDim*m_data.m_nbrEigenvalues*sizeof(T));
                        }
                    }
                    f_out.close();
                }
            }
            if(0 != eigenvectorsBuffer) 
            {
                delete[] eigenvectorsBuffer;
            }
            return;
        }

        //!
        //! Import all stored eigenvalue and eigenvector data from a text file
        //! and distribute data over program nodes
        //!
        void EigensystemFromFile(
            const std::string fileName,     //!<    Name of input file
            bool readValues,                //!<    Optionally read eigenvalues
            bool readVectors,               //!<    Optionally read eigenvectors
            const io::fileFormat_t format,  //!<    Format of file
            utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        {
            double* eigenvaluesBuffer = 0;
            T* eigenvectorsBuffer = 0;
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                //  Set flag to include eigenvalues with eigenvectors
                if(readVectors)     
                {
                    readValues = readVectors;
                }
                utilities::cout.AdditionalInfo()<<"\n\t- READING EIGENSYSTEM DATA FROM FILE "<<fileName.c_str()<<std::endl;
                std::ifstream f_in;
                utilities::GenFileStream(f_in, fileName, format, mpi);
                if(!mpi.m_exitFlag)
                {
                    if(readValues)
                    {
                        if(io::_TEXT_ == format)
                        {
                            std::string line;
                            while(getline(f_in, line))
                            {
                                if(line=="##  NBR EIGENVALUES") break;
                            }
                            f_in >> m_data.m_nbrEigenvalues;
                            eigenvaluesBuffer = new (std::nothrow) double[m_data.m_nbrEigenvalues];
                            while(getline(f_in, line))
                            {
                                if(line=="##  EIGENVALUES")	break;
                            }
                            for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                            {
                                f_in >> eigenvaluesBuffer[i];
                            }
                        }
                        else if(io::_BINARY_ == format)
                        {
                            f_in.read(reinterpret_cast<char*>(&(m_data.m_nbrEigenvalues)), sizeof(unsigned int));
                            eigenvaluesBuffer = new (std::nothrow) double[m_data.m_nbrEigenvalues];
                            f_in.read(reinterpret_cast<char*>(eigenvaluesBuffer), m_data.m_nbrEigenvalues*sizeof(double));
                        }
                        m_data.m_eigenvaluesCalculated = true;
                    }
                    if(readVectors)
                    {
                        if(io::_TEXT_ == format)
                        {
                            std::string line;
                            while(getline(f_in, line))
                            {
                                if(line=="##  FOCK SPACE DIMENSION") break;
                            }       
                            f_in >> m_data.m_fockSpaceDim;
                            eigenvectorsBuffer = new (std::nothrow) T[m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim];
                            while(getline(f_in,line))
                            {
                                if(line=="##  EIGENVECTORS") break;
                            }
                            T* p_eigs = eigenvectorsBuffer;
                            for(iSize_t i=0; i<m_data.m_nbrEigenvalues; ++i)
                            {
                                for(fock_t j=0; j<m_data.m_fockSpaceDim; ++j, ++p_eigs)
                                {
                                    utilities::FromStream(f_in, *p_eigs);
                                }           
                            }
                        }
                        else if(io::_BINARY_ == format)
                        {
                            f_in.read(reinterpret_cast<char*>(&(m_data.m_fockSpaceDim)), sizeof(unsigned int));
                            eigenvectorsBuffer = new (std::nothrow) T[m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim];
                            f_in.read(reinterpret_cast<char*>(eigenvectorsBuffer), m_data.m_nbrEigenvalues*m_data.m_fockSpaceDim*sizeof(T));
                        }
                    }
                    f_in.close();
                    m_data.m_eigenvectorsCalculated = true;
                }
            }
            //  MPI sync parameters required on all nodes
            mpi.Sync(&m_data.m_eigenvaluesCalculated, 1, 0);
            mpi.Sync(&m_data.m_eigenvectorsCalculated, 1, 0);
            mpi.Sync(&m_data.m_nbrEigenvalues, 1, 0);
            if(m_data.m_eigenvaluesCalculated)
            {
                m_eigenvalues.resize(m_data.m_nbrEigenvalues);
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    memcpy(m_eigenvalues.data(), eigenvaluesBuffer, m_data.m_nbrEigenvalues*sizeof(double));
                }
                MPI_Barrier(mpi.m_comm);
                //  Set the same eigenvalues on each node
                mpi.Sync(&m_eigenvalues, 0);
                //  Clear the buffer
                if(0 != eigenvaluesBuffer)  
                {
                    delete[] eigenvaluesBuffer;
                }
            }
            if(m_data.m_eigenvectorsCalculated)
            {
                //  Call MPI scatter on the eigenvectors
                mpi.DivideTasks(mpi.m_id,m_data.m_fockSpaceDim, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false); 
                fock_t eigenvectorNodeDim = mpi.m_lastTask - mpi.m_firstTask + 1;
                m_eigenvectors.resize(eigenvectorNodeDim*m_data.m_nbrEigenvalues);
                for(iSize_t e=0; e<m_data.m_nbrEigenvalues; ++e)
                {
                    MPI_Status status;
                    mpi.Scatter<T>(eigenvectorsBuffer+e*m_data.m_fockSpaceDim, m_data.m_fockSpaceDim,
                                   m_eigenvectors.data()+e*eigenvectorNodeDim, eigenvectorNodeDim, 0, mpi.m_comm, status);
                }
                if(0!=eigenvectorsBuffer)
                {
                    delete[] eigenvectorsBuffer;
                }
            }
            return;
        }
    };
}   //  End namespace diagonalization
#endif
