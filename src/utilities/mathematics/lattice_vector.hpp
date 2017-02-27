////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a class template to represent 2D lattice vectors. 
//!     It also implements a class template to contain a list of lattice vectors
//!     satisfying momentum conservation up to specified reciprocal lattice
//!     vectors.
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

#ifndef _LATTICE_VECTOR_HPP_INCLUDED_
#define _LATTICE_VECTOR_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../general/template_tools.hpp"
#include "../wrappers/mpi_wrapper.hpp"
#include "modulo.hpp"
#include <vector>

namespace utilities
{
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief Define a data type to contain a 2D vector, used to represent 
    //! lattice vectors. This data type has various arithmetic operator overloads 
    //! declared.
    //!
    //////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    class LatticeVector2D
    {
        private:
        int m_yDim;    //!<    Leading y array dimension for a 2D array
        T m_kx;        //!<    x-component of lattice vector
        T m_ky;        //!<    y-component of lattice vector
        public:
        //!
        //! Default constructor - initialize all members to zero
        //!
        LatticeVector2D()
        :
            m_yDim(0),
            m_kx(0),
            m_ky(0)
        {}
        //!
        //! A constructor to set the leading array dimension and initialize
        //! ky and kx to zero
        //!
        LatticeVector2D(
            const int yDim)    //!<    leading dimension of array
        :
            m_kx(0),
            m_ky(0)
        {
            m_yDim = yDim;
        }
        //!
        //! Constructor from a pair of values
        //!
        LatticeVector2D(
            const T kx,    //!<    set kx value
            const T ky)    //!<    set ky value
        :
            m_yDim(0),
            m_kx(kx),
            m_ky(ky)
        {}
        
        //!
        //! Copy constructor implementation
        //!
        LatticeVector2D(
            const LatticeVector2D& other)    //!<    Another instance of the object
        :
            m_yDim(other.m_yDim),
            m_kx(other.m_kx),
            m_ky(other.m_ky)
        {}

        //!
        //! Specialised copy constructor with a rescaling factor
        //!
        LatticeVector2D(
            const LatticeVector2D& other,   //!<    Another instance o fthe object
            const T scaleFactor)            //!<    rescaling factor
        :
            m_yDim(other.m_yDim),
            m_kx(other.m_kx*scaleFactor),
            m_ky(other.m_ky*scaleFactor)
        {}
               
        //!
        //! Overload of the equals (=) operator
        //!
        LatticeVector2D& operator=(const LatticeVector2D& other)
        {
            m_yDim = other.m_yDim;
            m_kx = other.m_kx;
            m_ky = other.m_ky;
            return *this;
        }
     
        //!
        //! Overload of the accumulation (+=) operator
        //!
        LatticeVector2D& operator+=(const LatticeVector2D& other)
        {
            if(m_yDim != other.m_yDim)
            {
                std::cerr<<"ERROR WITH LinearMomentum object - leading array dimension mismatch!"<<std::endl;
                exit(EXIT_FAILURE);
            }
            m_kx += other.m_kx;
            m_ky += other.m_ky;
            return *this;
        }
        
        //!
        //! Overload of the accumulation (-=) operator
        //!
        LatticeVector2D& operator-=(const LatticeVector2D& other)
        {
            if(m_yDim != other.m_yDim)
            {
                std::cerr<<"ERROR WITH LinearMomentum object - leading array dimension mismatch!"<<std::endl;
                exit(EXIT_FAILURE);
            }
            m_kx -= other.m_kx;
            m_ky -= other.m_ky;
            return *this;
        }
          
        //!
        //! Define a modulo operator to act on both kx and ky elements in the class 
        //!
        void Modulo(
            const T kxMod,    //!<    Modulo for x component
            const T kyMod)    //!<    Modulo for y component
        {
            m_kx = utilities::Modulo(m_kx, kxMod);
            m_ky = utilities::Modulo(m_ky, kyMod);
            m_kx = m_kx<0 ? m_kx + kxMod : m_kx;
            m_ky = m_ky<0 ? m_ky + kyMod : m_ky;
            return;
        }
         
        //!
        //! Define an equals operator to compare with a given pair of kx,ky values
        //!
        bool Equals(
            const T kx,       //!<    Comparison value for x component
            const T ky)       //!<    Comparison value for y component
            const
        {
            return (m_kx == kx) && (m_ky == ky);
        }
                  
        //!
        //! Map an array index i to a pair of fractional reciprocal lattice 
        //! indices x,y represented as co-ordiantes on a 2D k-space grid by 
        //! defining an appropriate constructor 
        //!
        void SetFromOrbital(
            const int i,      //!<    combined orbital index
            const int yDim)   //!<    leading dimension of array
        {
            m_yDim = yDim;
            m_ky = utilities::Modulo(i, m_yDim);
            m_kx = floor((double)i/(m_yDim));
        }
            
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Map a set of coordinates on a 2D k-space grid back to a single 
        //! orbital index
        //!
        //! This function performs the exact inverse operation of SetFromOrbital
        //!
        ////////////////////////////////////////////////////////////////////////////////
        int GetOrbital() const
        {
            return (int) (m_kx*m_yDim+m_ky);
        }
        //!
        //! Return kx parameter
        //!
        T GetKx() const
        {
            return m_kx;
        }       
        //!
        //! Return ky parameter
        //!
        T GetKy() const
        {
            return m_ky;
        }
        //!
        //! Return kx parameter
        //!
        void SetKx(const T kx)
        {
            m_kx = kx;
        }
        //!
        //! Return ky parameter
        //!
        void SetKy(const T ky)
        {
            m_ky = ky;
        }               
        //!
        //! Synchronize the class when run in parallel
        //!
        void MpiSync(
            const int nodeId,                   //!<    Node to sync with
            const utilities::MpiWrapper& mpi)   //!<    Instance of mpi wrapper class
        {
            MPI_Bcast(&m_kx, 1, mpi.GetType<T>(), nodeId, mpi.m_comm);
            MPI_Bcast(&m_ky, 1, mpi.GetType<T>(), nodeId, mpi.m_comm);
            MPI_Bcast(&m_yDim, 1, mpi.GetType<int>(), nodeId, mpi.m_comm);
        } 
        //!
        //! Print out the values contained
        //!
        void Print() const
        {
            std::cout<<"\tkx= "<<(int)m_kx<<"\tky="<<(int)m_ky<<std::endl;
        }
    };
};  //  End namespace utilities
#endif
