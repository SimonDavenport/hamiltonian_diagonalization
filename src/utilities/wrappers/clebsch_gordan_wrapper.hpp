////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a wrapper for calculation of Clebsch-Gordan 
//!     coefficients. Coefficients are calculated using the algorithm described 
//!     here: arXiv:1009.0437v1. See:
//!     http://homepages.physik.uni-muenchen.de/~vondelft/Papers/ClebschGordan/
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

#ifndef _CLEBSCH_GORDAN_WRAPPER_HPP_INCLUDED_
#define _CLEBSCH_GORDAN_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../mathematics/clebsch_gordan_coefficients.hpp"
#include "../data_structures/multi_key_hash.hpp"
#include "../general/cout_tools.hpp"
#if _DEBUG_
#include "../general/debug.hpp"
#endif
typedef uint8_t j_type;
typedef int16_t m_type;

namespace utilities
{
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief The ClebschGordonCoefficients class contains a look-up table of 
    //! pre-calculated Clebsch-Gordan coefficients.
    //!
    //! NOTE: All j and m angular momentum quantum numbers are multiplied by 2 
    //! by convention in order to store as integer rather than half integer values
    //!
    //! NOTE: The maximum allowed j value is 256 (which means the maximum quantum
    //! number value is 128).
    //////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    class ClebschGordanCoefficients
    {
        private:
        utilities::MultiHashMap<double> m_map;

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief  Convert signed int type used to represent the z-component of
        //! angular momentum index into an unsigned type suitable for bitwise operations
        ////////////////////////////////////////////////////////////////////////////////
        j_type ConvertToUnsigned(
            m_type value,       //!<    Integer value to be converted to an unsigned type
            j_type maxJ)        //!<    Maximum j value to shift by    
            const
        {
            return (value+maxJ)/2;
        }

        public:
        //!
        //! \brief  Default constructor declared empty 
        //!
        ClebschGordanCoefficients()
        {}
        //!
        //! \brief Destructor
        //!
        ~ClebschGordanCoefficients()
        {
            this->Clear();
        }
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief  Generate a set of tables of Clebsch-Gordan coefficients
        //! for multiplet labels j1 and j2 only and add to the current hash
        //! table of values. 
        //!
        //! NOTE: the maximum quanutm number value is 128 (represented by input j1=256).
        //! Any larger number will result in a truncation when converted to a j_type
        //! variable
        ////////////////////////////////////////////////////////////////////////////////
        void AddMultiplet(
            const j_type j1,       //!<    j1 between 0 and 256 represeting 2*j1
            const j_type j2)       //!<    j1 between 0 and 256 represeting 2*j2
        {
            //  Call the "clebsch" code to calculate all relevant SU(2) CG coefficients
            //  then reassign the CG coefficients to a more convenient hash map 
            //  indexes by conventional SU(2) weights and roots. 
            clebsch::Weight S(2);
            S(1) = j1;
            S(2) = 0;
            clebsch::Weight Sprime(2);
            Sprime(1) = j2;
            Sprime(2) = 0;
            clebsch::Decomposition decomp(S, Sprime);
            for(j_type d = 0; d < decomp.Size(); ++d) 
            {
                const clebsch::Coefficients C(decomp(d),S, Sprime);
                j_type dimS = S.Dimension();
                j_type dimSprime = Sprime.Dimension();
                j_type dimSdoubleprime = decomp(d).Dimension();
                for(j_type m = 0; m < C.multiplicity; ++m)
                {
                    for(j_type m3 = 0; m3 < dimSdoubleprime; ++m3) 
                    {
                        for(j_type m1 = 0; m1 < dimS;++m1) 
                        {
                            for(j_type m2 = 0; m2 < dimSprime; ++m2) 
                            {
                                double x = double(C(m1, m2, m, m3));
                                if(fabs(x) > clebsch::EPS)
                                {
                                    j_type j3 = decomp(d)(1);
                                    m_map.Insert(utilities::Key(j1, m1, j2, m2, j3, m3)) = x;
                                }
                            }
                        }
                    }
                }
            }
            return;
        }

        //!
        //! \brief  Clear all currently stored coefficients
        //! 
        void Clear()
        {
            m_map.Clear();
        }
        //!
        //! \brief  Get the size of the table of non-zero coefficients
        //! 
        unsigned int GetSize()
        {
            return m_map.Size();
        }
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief  Return the value of the specified Clebsch-Gordan coefficient from
        //! a pre-tabulated value
        //!
        //! NOTE: All j and m angular momentum quantum numbers are multiplied by 2 
        //! by convention in order to store as integer rather than half integer values
        //!
        //! NOTE1: This funciton does not check to see if the multiplet j1,j2 has been
        //! calculated or not. If it has not yet been calcualted then there are 
        //! no values so far in the look-up table
        ////////////////////////////////////////////////////////////////////////////////
        double Value(
            const j_type j1,    //!<    j1 between 0 and 256 represeting 2*j1
            const m_type m1,    //!<    m1 between -128 and 128 represeting 2*m1
            const j_type j2,    //!<    j2 between 0 and 256 represeting 2*j1
            const m_type m2,    //!<    m2 between -128 and 128 represeting 2*m2
            const j_type j3,    //!<    j1 between 0 and 256 represeting 2*j1
            const m_type m3)    //!<    m3 between -128 and 128 represeting 2*m3
            const
        {                      
            j_type unsignedM1 = this->ConvertToUnsigned(m1,j1);
            j_type unsignedM2 = this->ConvertToUnsigned(m2,j2);
            j_type unsignedM3 = this->ConvertToUnsigned(m3,j3);
            return m_map.Value(utilities::Key(j1,unsignedM1,j2,unsignedM2,j3,unsignedM3));
        }
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    };  //  End class ClebschGordanCoefficients
}       //  End namespace utilities
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
#endif
