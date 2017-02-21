////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                
//!                                                                             
//!	 \file
//!     This file defines a class to store data structures describing 
//!     matrix elements of a Hamiltonian	           
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
#include "lookup_tables.hpp"

namespace diagonalization
{
    //////  QuadraticLookUpTables IMPLEMENTATION        ////////////////////////

    //!
    //! Update the stored dimension of the quadratic term look-up table
    //!
    iSize_t QuadraticLookUpTables::CalculateDim(
        const kState_t kMax)    //!<    Number of k states
        const
    {
        return kMax;
    }

    //!
    //! Find the k1 value corresponding to a given k2.
    //!
    void QuadraticLookUpTables::GetK1(
        kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
        iSize_t& nbrK1,             //!<    Set the number of returned values 
        const kState_t k2)          //!<    k2 index
        const
    {
        nbrK1 = 1;
        kRetrieveBuffer[0] =  k2;
    }

    //!
    //! Return a quadratic term coefficient for a given k1,k2
    //!
    dcmplx QuadraticLookUpTables::GetEkk(
        const kState_t k1,         //!<    k1 index
        const kState_t k2)         //!<    k2 index
        const
    {
        return m_vTable[k1];
    }

    //////  QuarticLookUpTables IMPLEMENTATION        //////////////////////////

    //!
    //! Update the stored dimension of the quadratic term look-up table
    //!
    iSize_t QuarticLookUpTables::CalculateDim(
        const kState_t kMax)    //!<    Number of k states
        const
    {
        //  In the quartic term table we store V_{k1,k2,k3,k4}, but with one of
        //  the k indices completely fixed by quantum number conservation
        return (kMax+1)*std::pow(kMax,2)/2;
    }

    //!
    //! Get the k1 value corresponding to a given k2,k3,k4 (where each k 
    //! represents a k1,ky pair).
    //!
    void QuarticLookUpTables::GetK1(
        kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
        iSize_t& nbrK1,             //!<    Set the number of returned values 
        const kState_t k2,          //!<    k2 index
        const kState_t k3,          //!<    k3 index
        const kState_t k4)          //!<    k4 index
        const
    {
        kState_t temp1 = k4;
        kState_t temp2 = k3;
        if(k4>k3)
        {
            //  Swap and refer to the k3, k4 term in the table
            temp1 = k3;
            temp2 = k4;
        }
        nbrK1 = 1;
        kRetrieveBuffer[0] =  m_kTable[(temp1*(m_kMax) + temp2 - (temp1*(temp1+1))/2 )*m_kMax + k2];
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Return a quartic term coefficient specified by a given k1,k2,k3,k4
    //!
    //! NOTE: This function is ONLY valid when the k1,k2,k3,k4 indices form
    //! a momentum conserving set - otherwise the coefficient does not exist
    //! e.g. V_0001 does not correspond to a momentum conserving set of k1,k2,k3,k4
    //! , but the GetVkkkk function will return a non-zero result nevertheless
    //! , and that result will be incorrect!
    //!
    ////////////////////////////////////////////////////////////////////////////////
    dcmplx QuarticLookUpTables::GetVkkkk(
        const kState_t k1,         //!<    k1 index
        const kState_t k2,         //!<    k2 index
        const kState_t k3,         //!<    k3 index
        const kState_t k4)         //!<    k4 index
        const
    {
        kState_t temp1 = k4;
        kState_t temp2 = k3;
        kState_t temp3 = k2;
        if(k4>k3)
        {
            //  Swap and refer to the k3, k4 term in the table (also requires swapping k1 and k2)
            temp1 = k3;
            temp2 = k4;
            temp3 = k1;
        }
        return m_vTable[(temp1*(m_kMax) + temp2 - (temp1*(temp1+1))/2 )*m_kMax + temp3];
    }

    //////  QuadraticLookUpHashTables IMPLEMENTATION        ////////////////////

    //!
    //! Get the highest number of k1 values that will be returned
    //! for a given k2. This allows us to pre-allocate a memory buffer
    //! to store the return values statically
    //!
    
    iSize_t QuadraticLookUpHashTables::GetMaxKCount() const
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
    void QuadraticLookUpHashTables::GetK1(
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
    dcmplx QuadraticLookUpHashTables::GetEkk(
        const kState_t k1,  //!<    k1 value
        const kState_t k2)  //!<    k2 value
        const
    {
        return m_vTable.Value(utilities::Key(k1, k2));
    }

    //////  QuarticLookUpHashTables IMPLEMENTATION        //////////////////////

    //!
    //! Get the highest number of k1 values that will be returned
    //! for a given k2,k3,k4. This allows us to pre-allocate a memory buffer
    //! to store the return values statically
    //!
    iSize_t QuarticLookUpHashTables::GetMaxKCount() const
    {
        return m_kTable.GetMaxCount();
    }

    //!
    //! Get the k1 value for a given k2,k3,k4 with a 4 momentum conservation 
    //! law
    //!
    void QuarticLookUpHashTables::GetK1(
        kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
        iSize_t& nbrK1,             //!<    Set the number of returned values             
        const kState_t k2,          //!<    Given k2 value
        const kState_t k3,          //!<    Given k3 value
        const kState_t k4)          //!<    Given k4 value
        const
    {
        //  Get the corresponding list of k1 values
        m_kTable.Value(kRetrieveBuffer,nbrK1,utilities::Key(k2, k3, k4));
    }

    //!
    //! Return a specified coefficient of the quartic term
    //!
    dcmplx QuarticLookUpHashTables::GetVkkkk(
        const kState_t k1,  //!<    Given k1 value
        const kState_t k2,  //!<    Given k2 value
        const kState_t k3,  //!<    Given k3 value
        const kState_t k4)  //!<    Given k4 value
        const
    {
        return m_vTable.Value(utilities::Key(k1,k2,k3,k4)); 
    }
}   //  End namespace diagonalization
