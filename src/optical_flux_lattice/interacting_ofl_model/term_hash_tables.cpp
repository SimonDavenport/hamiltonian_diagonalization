////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store data structures describing 
//!     terms of a Hamiltonian	           
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
#include "term_hash_tables.hpp"

namespace diagonalization
{
    //////  QuadraticTermHashTables IMPLEMENTATION        ////////////////////

    //!
    //! Get the highest number of k1 values that will be returned
    //! for a given k2. This allows us to pre-allocate a memory buffer
    //! to store the return values statically
    //!
    iSize_t QuadraticTermHashTables::GetMaxKCount() const
    {
        iSize_t maxCount = 1;
        for(kState_t k2=0; k2<m_kMax; ++k2)
        {
            maxCount = std::max(maxCount, (iSize_t)m_kTable.Count(k2));
        }
        return maxCount;
    }

    //!
    //! Get the k1 value for a given k2, with a 2 value conservation law
    //!
    void QuadraticTermHashTables::GetK1(
        kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
        iSize_t& nbrK1,             //!<    Set the number of returned values             
        const kState_t k2)          //!<    Given k2 value
        const
    {
        //  Get the corresponding list of k1 values
        m_kTable.Value(kRetrieveBuffer, nbrK1, k2);
    }

    //!
    //! Return a specified coefficient of the quadratic term
    //!
    dcmplx QuadraticTermHashTables::GetEkk(
        const kState_t k1,  //!<    k1 value
        const kState_t k2)  //!<    k2 value
        const
    {
        return m_vTable.Value(utilities::Key(k1, k2));
    }
    
    //!
    //! Write quadratic terms to a file
    //!
    void QuadraticTermHashTables::ToFile(
        const std::string fileName, //!<    Name of file
        const std::string format,   //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
    {
        this->ToFileBase(fileName, format, 2, mpi);
    }
    
    //!
    //! Read quadratic terms from a file
    //!
    void QuadraticTermHashTables::FromFile(
        const std::string fileName, //!<    Name of file
        const std::string format,   //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
    {
        this->FromFileBase(fileName, format, mpi);
    }

    //////  QuarticLookUpHashTables IMPLEMENTATION        //////////////////////

    //!
    //! Get the highest number of k1 values that will be returned
    //! for a given k2,k3,k4. This allows us to pre-allocate a memory buffer
    //! to store the return values statically
    //!
    iSize_t QuarticTermHashTables::GetMaxKCount() const
    {
        return m_kTable.GetMaxCount();
    }

    //!
    //! Get the k1 value for a given k2,k3,k4 with a 4 momentum conservation 
    //! law
    //!
    void QuarticTermHashTables::GetK1(
        kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
        iSize_t& nbrK1,             //!<    Set the number of returned values             
        const kState_t k2,          //!<    Given k2 value
        const kState_t k3,          //!<    Given k3 value
        const kState_t k4)          //!<    Given k4 value
        const
    {
        m_kTable.Value(kRetrieveBuffer, nbrK1, utilities::Key(k2, k3, k4));
    }

    //!
    //! Return a specified coefficient of the quartic term
    //!
    dcmplx QuarticTermHashTables::GetVkkkk(
        const kState_t k1,  //!<    Given k1 value
        const kState_t k2,  //!<    Given k2 value
        const kState_t k3,  //!<    Given k3 value
        const kState_t k4)  //!<    Given k4 value
        const
    {
        return m_vTable.Value(utilities::Key(k1, k2, k3, k4)); 
    }
    
    //!
    //! Write quartic terms to a file
    //!
    void QuarticTermHashTables::ToFile(
        const std::string fileName, //!<    Name of file
        const std::string format,   //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
    {
        this->ToFileBase(fileName, format, 4, mpi);
    }
    
    //!
    //! Read quartic terms from a file
    //!
    void QuarticTermHashTables::FromFile(
        const std::string fileName, //!<    Name of file
        const std::string format,   //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
    {
        this->FromFileBase(fileName, format, mpi);
    }
}   //  End namespace diagonalization
