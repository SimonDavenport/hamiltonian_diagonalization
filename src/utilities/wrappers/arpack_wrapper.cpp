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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "arpack_wrapper.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

namespace utilities
{
    namespace linearAlgebra
    {
        //!
        //! Default constructor with initializer list for default parameters
        //!
        ArpackWrapper::ArpackWrapper()
        :   
            m_typeOfProblem('I'),
            m_eigenvalueMagnitude('S'),
            m_whatMagnitude('R'),
            m_upperOrLowerMatrix('U'),
            m_maxIterations(2000),
            m_provideInitial(0),
            m_storeFinalVector(false),
            m_initialVectorFile("initial_vector.bin"),
            m_finalVectorFile("initial_vector.bin")
        {}
        //!
        //! Copy constructor
        //!
        ArpackWrapper::ArpackWrapper(const ArpackWrapper& other)
        :   
            m_typeOfProblem(other.m_typeOfProblem),
            m_eigenvalueMagnitude(other.m_eigenvalueMagnitude),
            m_whatMagnitude(other.m_whatMagnitude),
            m_upperOrLowerMatrix(other.m_upperOrLowerMatrix),
            m_maxIterations(other.m_maxIterations),
            m_provideInitial(other.m_provideInitial),
            m_storeFinalVector(other.m_storeFinalVector),
            m_initialVectorFile(other.m_initialVectorFile),
            m_finalVectorFile(other.m_finalVectorFile)
        {}
        //!
        //! Assignment operator
        //!
        ArpackWrapper& ArpackWrapper::operator=(const ArpackWrapper& other)
        {
            m_typeOfProblem = other.m_typeOfProblem;
            m_eigenvalueMagnitude = other.m_eigenvalueMagnitude;
            m_whatMagnitude = other.m_whatMagnitude;
            m_upperOrLowerMatrix = other.m_upperOrLowerMatrix;
            m_maxIterations = other.m_maxIterations;
            m_provideInitial = other.m_provideInitial;
            m_storeFinalVector = other.m_storeFinalVector;
            m_initialVectorFile = other.m_initialVectorFile;
            m_finalVectorFile = other.m_finalVectorFile;
            return *this;
        }
        //!
        //! Synchronize parameter values with those on the master node
        //!
        void ArpackWrapper::MpiSync(
            const int nodeId,       //!<    Node to synchronize with
            const MpiWrapper& mpi)  //!<    Instance of the MPI wrapper class
        {
            mpi.Sync(&m_typeOfProblem, 1, nodeId);
            mpi.Sync(&m_eigenvalueMagnitude, 1, nodeId);
            mpi.Sync(&m_whatMagnitude, 1, nodeId);
            mpi.Sync(&m_upperOrLowerMatrix, 1, nodeId);
            mpi.Sync(&m_maxIterations, 1, nodeId);
            mpi.Sync(&m_provideInitial, 1, nodeId);
            mpi.Sync(&m_storeFinalVector, 1, nodeId);
            mpi.Sync(m_initialVectorFile, nodeId);
            mpi.Sync(m_finalVectorFile, nodeId);
        }
        
        //!
        //! Set the type of problem: 'I' for standard eigenvalue problem or
        //! 'G' for generalized eigenvalue problem
        //!
        void ArpackWrapper::SetProblemType(
            const char typeOfProblem)   //!<    Set type of problem
        {
            if('I' != typeOfProblem && 'G' != typeOfProblem)
            {
                std::cerr<<"ERROR in SetProblemType: argument "<<typeOfProblem;
                std::cerr<<" not recognized!"<<std::endl;
                return;
            }
            else
            {
                m_typeOfProblem = typeOfProblem;
            }
        }
        
        //!
        //! Identify the part of the spectrum to investigate: 'S' sets the 
        //! smallest eigenvalues, 'L' the largest.See also "SetWhatMagnitude" function
        //!
        void ArpackWrapper::SetEigenvalueMagnitude(
            const char eigenvalueMagnitude) //!<    Set eigenvalue magnitude
        {
            if('S' != eigenvalueMagnitude && 'L' != eigenvalueMagnitude)
            {
                std::cerr<<"ERROR in SetEigenvalueMagnitude: argument "<<eigenvalueMagnitude;
                std::cerr<<" not recognized!"<<std::endl;
                return;
            }
            else
            {
                m_eigenvalueMagnitude = eigenvalueMagnitude;
            }
        }
        
        //!
        //! Identify the part of the spectrum to investigate: 'I' means we look
        //! at eigenvalues with the 'S' or 'L' imaginary part, 'R' for real part
        //! and 'M' for the absolute value. See also "SetEigenvalueMagnitude" function
        //!
        void ArpackWrapper::SetWhatMagnitude(
            const char whatMagnitude)       //!<    Set what the magnitude refers to 
        {
            if('I' != whatMagnitude && 'R' != whatMagnitude && 'M' != whatMagnitude)
            {
                std::cerr<<"ERROR in SetWhatMagnitude: argument "<<whatMagnitude;
                std::cerr<<" not recognized!"<<std::endl;
                return;
            }
            else
            {
                m_whatMagnitude = whatMagnitude;
            }
        }
        
        //!
        //! Specify whether the matrix passed to the ARPACK routine is stored in
        //! upper ('U') or lower ('L') triangular format. 
        //!
        void ArpackWrapper::SetUpperOrLower(
            const char upperOrLower)        //!<    Upper or lower triangular?
        {
            if('U' != upperOrLower && 'L' != upperOrLower)
            {
                std::cerr<<"ERROR in SetUpperOrLower: argument "<<upperOrLower;
                std::cerr<<" not recognized!"<<std::endl;
                return;
            }
            else
            {
                m_upperOrLowerMatrix = upperOrLower;
            }
        }
        
        //!
        //! Return the currently set value of the m_upperOrLowerMatrix parameter
        //!
        char ArpackWrapper::GetUpperOrLower()
            const
        {
           return m_upperOrLowerMatrix;
        }
        
        //!
        //! Set the maximum number of Lanczos iterations to perform before 
        //! quitting
        //!
        void ArpackWrapper::SetMaxIterations(
            const int maxIterations)        //!<    Set maximum iterations
        {
            m_maxIterations = maxIterations;
        }
        
        //!
        //! Set option to use initial Lanczos vector by specifying a file where
        //! it is stored. ARPACK allows one to define an initial vector, but by default
        //! starts with a random vector
        //!
        void ArpackWrapper::UseInitialVector()
        {
            m_provideInitial = 1;
        }
        
        //!
        //! \brief Set the initial vector file name
        //!
        void ArpackWrapper::SetInitialVectorFile(
            const std::string initialVectorFile)    //!<    Name of initial vector file
        {
            m_initialVectorFile = initialVectorFile;
        }
        
        //!
        //! Set option to store the first of the final Lanczos vectors, 
        //! potentially to be used as a new initial vector
        //!
        void ArpackWrapper::StoreFinalVector()
        {
            m_storeFinalVector = true;
        }
        
        //!
        //! \brief Set the final vector file name
        //!
        void ArpackWrapper::SetFinalVectorFile(
            const std::string finalVectorFile)  //!<    Name of final vector file
        {
             m_finalVectorFile = finalVectorFile;
        }

    }   //  End namespace linearAlgebra 
}   //  End namespace utilities 
