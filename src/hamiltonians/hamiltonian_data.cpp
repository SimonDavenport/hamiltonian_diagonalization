////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 03/02/2015
//!
//!  \file
//!		This file contains an underlying data structure to store variables in the
//!     any Hamiltonian class
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "hamiltonian_data.hpp"

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Write a value to a file for double types
//!
////////////////////////////////////////////////////////////////////////////////

void ValueToFile(std::ofstream& f_out,double* value)
{
    f_out.precision(15);
    f_out<<*value<<"\n";
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Write a value to a file for dcmplx types
//!
////////////////////////////////////////////////////////////////////////////////

void ValueToFile(std::ofstream& f_out,dcmplx* value)
{
    f_out.precision(15);                
    f_out<<value->real()<<" ";
    
    f_out.precision(15);
    f_out<<value->imag()<<"\n";
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Read a value from a file for double types
//!
////////////////////////////////////////////////////////////////////////////////

void ValueFromFile(std::ifstream& f_in,double* value)
{
    f_in >> *value;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Read a value from a file for dcmplx types
//!
////////////////////////////////////////////////////////////////////////////////

void ValueFromFile(std::ifstream& f_in,dcmplx* value)
{
    double tmpReal;
    double tmpImag;

    f_in >> tmpReal;
    f_in >> tmpImag;
    
    *value = dcmplx(tmpReal,tmpImag);
}
                                                                             
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Default constructor declaration
    //!
    ////////////////////////////////////////////////////////////////////////////////

    HamiltonianData::HamiltonianData()
    : 
    m_nbrParticles(1),
    m_nbrOrbitals(1),
    m_highestOrbital(0),
    m_storageMethod(utilities::_SPARSE_MAPPED_),
    m_matrixAllocated(false),
    m_eigenvaluesCalculated(false),
    m_eigenvectorsCalculated(false),
    m_fockSpaceDim(0),
    m_nodeDim(0),
    m_nbrEigenvalues(0)
    {
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor with basic arguments
    //!
    ////////////////////////////////////////////////////////////////////////////////

    HamiltonianData::HamiltonianData(
        const iSize_t nbrParticles,                     //!<    Number of particles
        const iSize_t nbrOrbitals,                      //!<    Number of orbitals
        const utilities::storageMethod_t storageMethod) //!<    Matrix storage method
    : 
    m_nbrParticles(nbrParticles),
    m_nbrOrbitals(nbrOrbitals),
    m_highestOrbital(0),
    m_storageMethod(storageMethod),
    m_matrixAllocated(false),
    m_eigenvaluesCalculated(false),
    m_eigenvectorsCalculated(false),
    m_fockSpaceDim(0),
    m_nodeDim(0),
    m_nbrEigenvalues(0)
    {
        m_highestOrbital = utilities::binary::number1 << (m_nbrOrbitals-1);
    }
       
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Copy constructor declaration
    //!
    ////////////////////////////////////////////////////////////////////////////////

    HamiltonianData::HamiltonianData(const HamiltonianData& other)
    : 
    m_nbrParticles(other.m_nbrParticles),
    m_nbrOrbitals(other.m_nbrOrbitals),
    m_highestOrbital(other.m_highestOrbital),
    m_storageMethod(other.m_storageMethod),
    m_matrixAllocated(other.m_matrixAllocated),
    m_eigenvaluesCalculated(other.m_eigenvaluesCalculated),
    m_eigenvectorsCalculated(other.m_eigenvectorsCalculated),
    m_fockSpaceDim(other.m_fockSpaceDim),
    m_nodeDim(other.m_nodeDim),
    m_nbrEigenvalues(other.m_nbrEigenvalues)
    {}
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Copy assignment operator declaration
    //!
    ////////////////////////////////////////////////////////////////////////////////

    HamiltonianData& HamiltonianData::operator=(const HamiltonianData& other)
    {
        m_nbrParticles = other.m_nbrParticles;
        m_nbrOrbitals = other.m_nbrOrbitals;
        m_highestOrbital = other.m_highestOrbital;
        m_storageMethod = other.m_storageMethod;
        m_matrixAllocated = other.m_matrixAllocated;
        m_eigenvaluesCalculated = other.m_eigenvaluesCalculated;
        m_eigenvectorsCalculated = other.m_eigenvectorsCalculated;
        m_fockSpaceDim = other.m_fockSpaceDim;
        m_nodeDim = other.m_nodeDim;
        m_nbrEigenvalues = other.m_nbrEigenvalues;

        return *this;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Destructor
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    HamiltonianData::~HamiltonianData()
    {
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Reset class flags
    //!
    ////////////////////////////////////////////////////////////////////////////////
     
    void HamiltonianData::Clear()
    {          
        m_matrixAllocated           = false;
        m_eigenvaluesCalculated     = false;
        m_eigenvectorsCalculated    = false;
        m_fockSpaceDim              = 0;
        m_nodeDim                   = 0;
    }       

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Set Fock Space dimension and dimension if stored on multiple nodes
    //!
    //! Eigenvector node dimension is always distributed over all nodes
    //!
    ////////////////////////////////////////////////////////////////////////////////
    
    void HamiltonianData::SetDataDimensions(
        const fock_t fockSpaceDimension,    //!<    Total Fock space dimension
        const iSize_t nbrNodes,             //!<    Specify the number of nodes
                                            //!     that the Fock space will be
                                            //!     distributed over
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        m_fockSpaceDim = fockSpaceDimension;
        
        mpi.DivideTasks(mpi.m_id,m_fockSpaceDim,nbrNodes,&mpi.m_firstTask,&mpi.m_lastTask,false);
        
        if(mpi.m_id < (int)nbrNodes)
        {
            m_nodeDim = mpi.m_lastTask - mpi.m_firstTask + 1;
        }
        else
        {
            m_nodeDim = 0;
        }
        
        //mpi.DivideTasks(mpi.m_id,m_fockSpaceDim,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);

        utilities::cout.DebuggingInfo()<<"\n\tON NODE "<<mpi.m_id<<" SET MATRIX NODE DIMENSION TO "<<m_nodeDim<<std::endl;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End diagonalization namespace
  
