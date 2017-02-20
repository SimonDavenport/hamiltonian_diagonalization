////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!                         \date Last Modified: 10/01/2015
//!
//!  \file
//!     This file provides wrappers for several ARPACK and PARPACK routines.
//!     See here http://www.caam.rice.edu/software/ARPACK/ for documentation.
//!
//!                    Copyright (C) 2015 Simon C Davenport
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

//////////////////////////////////////////////////////////////////////////////// 
//! \brief A namespace to contain any functions and utilities that I have 
//! written for use with any c++ program.
//!
////////////////////////////////////////////////////////////////////////////////

namespace utilities
{
    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief A function Namespace for linear algebra routines
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    namespace linearAlgebra
    {

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Default constructor with initializer list for default parameters
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        ArpackWrapper::ArpackWrapper()
        :   
            m_typeOfProblem('I'),
            m_eigenvalueMagnitude('S'),
            m_whatMagnitude('R'),
            m_upperOrLowerMatrix('U'),
            m_shift(0.0),
            m_maxIterations(2000),
            m_mode(arpack::_STANDARD_),
            m_provideInitial(0),
            m_storeFinalVector(false),
            m_initialVectorFile("initial_vector.bin"),
            m_finalVectorFile("initial_vector.bin"),
            m_matrixPower(1)
        {
            m_coefficientA.resize(m_maxMatrixPower);
            m_coefficientB.resize(m_maxMatrixPower);
        
            //  Set default coefficient values to solve standard eigenvalue problem
            
            m_coefficientA[0] = 1.0;
            m_coefficientB[0] = 0.0;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Copy constructor
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        ArpackWrapper::ArpackWrapper(const ArpackWrapper& other)
        :   
            m_typeOfProblem(other.m_typeOfProblem),
            m_eigenvalueMagnitude(other.m_eigenvalueMagnitude),
            m_whatMagnitude(other.m_whatMagnitude),
            m_upperOrLowerMatrix(other.m_upperOrLowerMatrix),
            m_shift(other.m_shift),
            m_maxIterations(other.m_maxIterations),
            m_mode(other.m_mode),
            m_provideInitial(other.m_provideInitial),
            m_storeFinalVector(other.m_storeFinalVector),
            m_initialVectorFile(other.m_initialVectorFile),
            m_finalVectorFile(other.m_finalVectorFile),
            m_matrixPower(other.m_matrixPower),
            m_coefficientA(other.m_coefficientA),
            m_coefficientB(other.m_coefficientB)
        {}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Assignment operator
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        ArpackWrapper& ArpackWrapper::operator=(const ArpackWrapper& other)
        {
            m_typeOfProblem = other.m_typeOfProblem;
            m_eigenvalueMagnitude = other.m_eigenvalueMagnitude;
            m_whatMagnitude = other.m_whatMagnitude;
            m_upperOrLowerMatrix = other.m_upperOrLowerMatrix;
            m_shift = other.m_shift;
            m_maxIterations = other.m_maxIterations;
            m_mode = other.m_mode;
            m_provideInitial = other.m_provideInitial;
            m_storeFinalVector = other.m_storeFinalVector;
            m_initialVectorFile = other.m_initialVectorFile;
            m_finalVectorFile = other.m_finalVectorFile;
            m_matrixPower = other.m_matrixPower;
            m_coefficientA = other.m_coefficientA;
            m_coefficientB = other.m_coefficientB;
            
            return *this;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Synchronize parameter values with those on the master node
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void ArpackWrapper::MpiSync(
            const int nodeId,       //!<    Node to synchronize with
            const MpiWrapper& mpi)  //!<    Instance of the MPI wrapper class
        {
            mpi.Sync(&m_typeOfProblem,1,nodeId);
            mpi.Sync(&m_eigenvalueMagnitude,1,nodeId);
            mpi.Sync(&m_whatMagnitude,1,nodeId);
            mpi.Sync(&m_upperOrLowerMatrix,1,nodeId);
            mpi.Sync(&m_shift,1,nodeId);
            mpi.Sync(&m_maxIterations,1,nodeId);
            mpi.Sync(&m_provideInitial,1,nodeId);
            mpi.Sync(&m_storeFinalVector,1,nodeId);
            mpi.Sync(m_initialVectorFile,nodeId);
            mpi.Sync(m_finalVectorFile,nodeId);
            mpi.Sync(&m_matrixPower,1,nodeId);
            mpi.Sync(&m_coefficientA[0],m_coefficientA.size(),nodeId);
            mpi.Sync(&m_coefficientB[0],m_coefficientB.size(),nodeId);
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
      
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set the ARPACK mode (_STANDARD_ or _SHIFT_INVERT_)
        //!
        //////////////////////////////////////////////////////////////////////////////////
      
        void ArpackWrapper::SetMode(
            const mode_t mode)      //!<    Set the new mode
        {
            m_mode = mode;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set the value of the shift (used in shift-invert mode only)
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void ArpackWrapper::SetShift(
            const double shift) //!<    Set the new shift
        {
            m_shift = shift;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set the type of problem: 'I' for standard eigenvalue problem or
        //! 'G' for generalized eigenvalue problem
        //!
        //////////////////////////////////////////////////////////////////////////////////
    
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
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Identify the part of the spectrum to investigate: 'S' sets the 
        //! smallest eigenvalues, 'L' the largest.See also "SetWhatMagnitude" function
        //!
        //////////////////////////////////////////////////////////////////////////////////

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

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Identify the part of the spectrum to investigate: 'I' means we look
        //! at eigenvalues with the 'S' or 'L' imaginary part, 'R' for real part
        //! and 'M' for the absolute value. See also "SetEigenvalueMagnitude" function
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
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

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Specify whether the matrix passed to the ARPACK routine is stored in
        //! upper ('U') or lower ('L') triangular format. 
        //!
        //////////////////////////////////////////////////////////////////////////////////

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

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Return the currently set value of the m_upperOrLowerMatrix parameter
        //!
        //////////////////////////////////////////////////////////////////////////////////

        char ArpackWrapper::GetUpperOrLower()
            const
        {
           return m_upperOrLowerMatrix;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set the maximum number of Lanczos iterations to perform before 
        //! quitting
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void ArpackWrapper::SetMaxIterations(
            const int maxIterations)        //!<    Set maximum iterations
        {
            m_maxIterations = maxIterations;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set option to use initial Lanczos vector by specifying a file where
        //! it is stored. ARPACK allows one to define an initial vector, but by default
        //! starts with a random vector
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void ArpackWrapper::UseInitialVector()
        {
            m_provideInitial = 1;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set the initial vector file name
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void ArpackWrapper::SetInitialVectorFile(
            const std::string initialVectorFile)    //!<    Name of initial vector file
        {
            //  Set the file name of the initial vector file
            
            m_initialVectorFile = initialVectorFile;
        }
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set option to store the first of the final Lanczos vectors, 
        //! potentially to be used as a new initial vector
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void ArpackWrapper::StoreFinalVector()
        {
            //  Update the final vector flag
            
            m_storeFinalVector = true;
        }
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Set the final vector file name
        //!
        //////////////////////////////////////////////////////////////////////////////////
        
        void ArpackWrapper::SetFinalVectorFile(
            const std::string finalVectorFile)  //!<    Name of final vector file
        {
            //  Set the file name of the final vector file
            
             m_finalVectorFile = finalVectorFile;
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Specify the form of a multiplicative preconditioner. The general form
        //! is set to solve for the eigenvalues of (a1*A+b1)(a2*A+B2)...
        //! 
        //! The idea is to transform the eigenvalues spectrum to a form where those
        //! eigenvalues we want to calculated are well separated from the remainder. 
        //!
        //! This is done at the expense of having to perform more than one matrix-vector
        //! multiplication at each step. 
        //!
        //! One approach that often works when finding the lowest lying eigenvalues 
        //! is to determine the eigenvalues of -A^2.
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void ArpackWrapper::SetPreconditioner(
            const int matrixPower,      //!<    Set power of the matrix used
            double* a,                  //!<    List of preconditioner a coefficients
            double* b)                  //!<    List of preconditioner b coefficients
        {
            if(matrixPower>m_maxMatrixPower)
            {
                std::cerr<<"\n\tWARNING: matrix power set in ARPACK preconditioner > "<<m_maxMatrixPower<<" not allowed"<<std::endl;
                
                m_matrixPower = m_maxMatrixPower;
            }
            else
            {
                m_matrixPower = matrixPower;
            }

            for(int i=0;i<matrixPower;i++)
            {
                m_coefficientA[i] = a[i];
                m_coefficientB[i] = b[i];
            }
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
        //////////////////////////////////////////////////////////////////////////////////
        //! \brief Get the current values of the preconditioner parameters:
        //! m_matrixPower, m_coefficientA and m_coefficientB
        //!
        //////////////////////////////////////////////////////////////////////////////////

        void ArpackWrapper::GetPreconditioner(
            int& matrixPower,           //!<    Power of the matrix used
            std::vector<double>& a,     //!<    List of preconditioner a coefficients
            std::vector<double>& b)     //!<    List of preconditioner b coefficients
            const
        {
            matrixPower = m_matrixPower;
            a = m_coefficientA;
            b = m_coefficientB;    
        }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
    }   //  End namespace linearAlgebra 
    
}   //  End namespace utilities 
