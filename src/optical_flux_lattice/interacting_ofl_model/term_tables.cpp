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
#include "term_tables.hpp"

namespace diagonalization
{
    //////  QuadraticTermTables IMPLEMENTATION        ////////////////////////

    //!
    //! Update the stored dimension of the quadratic term look-up table
    //!
    iSize_t QuadraticTermTables::CalculateDim(
        const kState_t kMax)    //!<    Number of k states
        const
    {
        return kMax;
    }

    //!
    //! Find the k1 value corresponding to a given k2.
    //!
    void QuadraticTermTables::GetK1(
        kState_t* kRetrieveBuffer,  //!<    Buffer to store retrieved k values
        iSize_t& nbrK1,             //!<    Set the number of returned values 
        const kState_t k2)          //!<    k2 index
        const
    {
        nbrK1 = 1;
        kRetrieveBuffer[0] = k2;
    }

    //!
    //! Return a quadratic term coefficient for a given k1,k2
    //!
    dcmplx QuadraticTermTables::GetEkk(
        const kState_t k1,         //!<    k1 index
        const kState_t k2)         //!<    k2 index
        const
    {
        return m_vTable[k1];
    }

    //!
    //! Write quadratic terms to a file
    //!
    void QuadraticTermTables::ToFile(
        const std::string fileName, //!<    Name of file
        const std::string format,   //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
        const
    {
        this->ToFileBase(fileName, format, 2, mpi);
    }
    
    //!
    //! Read quadratic terms from a file
    //!
    void QuadraticTermTables::FromFile(
        const std::string fileName, //!<    Name of file
        const std::string format,   //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
    {
        this->FromFileBase(fileName, format, mpi);
    }

    //////  QuarticTermTables IMPLEMENTATION        //////////////////////////

    //!
    //! Update the stored dimension of the quadratic term look-up table
    //!
    iSize_t QuarticTermTables::CalculateDim(
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
    void QuarticTermTables::GetK1(
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

    //!
    //! Set the k1 value corresponding to a given k2, k3, k4
    //!
    void QuarticTermTables::SetK1(
        const kState_t k1,          //!<    k1 index to be set
        const kState_t k2,          //!<    k2 index
        const kState_t k3,          //!<    k3 index
        const kState_t k4)          //!<    k4 index
    {
        kState_t temp1 = k4;
        kState_t temp2 = k3;
        if(k4>k3)
        {
            //  Swap and refer to the k3, k4 term in the table
            temp1 = k3;
            temp2 = k4;
        }
        m_kTable[(temp1*(m_kMax) + temp2 - (temp1*(temp1+1))/2 )*m_kMax + k2] = k1;
    }  

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Return a quartic term coefficient specified by a given k1,k2,k3,k4
    //!
    //! NOTE: This function is ONLY valid when the k1,k2,k3,k4 indices form
    //! a momentum conserving set - otherwise the coefficient does not exist
    //! e.g. V_0001 does not correspond to a momentum conserving set of k1,k2,k3,k4
    //! , but the GetVkkkk function will return a non-zero result nevertheless
    //! , and that result will be incorrect!
    ////////////////////////////////////////////////////////////////////////////////
    dcmplx QuarticTermTables::GetVkkkk(
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
    
    //!
    //! Set the value of the term Vkkkk for a given k1, k2, k3, k4
    //!
    void QuarticTermTables::SetVkkkk(
        const dcmplx Vkkkk,        //!<    New Vkkkk value to set
        const kState_t k1,         //!<    k1 index
        const kState_t k2,         //!<    k2 index
        const kState_t k3,         //!<    k3 index
        const kState_t k4)         //!<    k4 index
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
        m_vTable[(temp1*(m_kMax) + temp2 - (temp1*(temp1+1))/2 )*m_kMax + temp3] = Vkkkk;
    }
    
    //!
    //! Write quartic terms to a file
    //!
    void QuarticTermTables::ToFile(
        const std::string fileName, //!<    Name of file
        std::string format,         //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
        const
    {
        this->ToFileBase(fileName, format, 4, mpi);
    }
    
    //!
    //! Read quartic terms from a file
    //!
    void QuarticTermTables::FromFile(
        const std::string fileName, //!<    Name of file
        std::string format,         //!<    Format of file (e.g. "binary", "text")
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class 
    {
        this->FromFileBase(fileName, format, mpi);
    }
}   //  End namespace diagonalization
