////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 27/03/2015              
//!                                                              
//!	 \file
//!     This file generates a set of two-body pseudopotential coefficients
//!     for the Coulomb interaction, implemented in the sphere geometry  
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
#include "coulomb_pseudopotentials.hpp"

namespace diagonalization
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Funciton to generate the set of two-body Haldane Pseudopotentials
    //! for the Coulomb interaciton in the Sphere goemetry. The analytic form
    //! for the lowest Landau level is
    //!
    //! V_L = 2/R [4Q-2L]C[2Q-L] * [4Q+2L+2]C[2Q+L+1] / ([4Q+2]C[2Q+1])^2
    //!
    //! Where R is the sphere raduis, given by sqrt(Q)
    //!
    ////////////////////////////////////////////////////////////////////////////////

    std::vector<double> GenCoulombPseudopotentials(
        const iSize_t nbrOrbitals,  //!<    Total number of orbitals
        const iSize_t llIndex,      //!<    Landau level index (0 is lowest)
        const double mulFactor)     //!<    Additional multiplicative factor
    {
        //  Calcualte 2*Maximum Lz value in the highest LL
        iSize_t maxLz = 0;
    
        if(llIndex==0)
        {
            maxLz = nbrOrbitals-1;
        }
        else
        {
            std::cerr<<"llIndex > 0 NOT YET PROGRAMMED!"<<std::endl;
        }
        //else if(llIndex==1)
        //{
        //    maxLz = (nbrOrbitals-2)/2+1;
        //}

        //   Vector containing pseudopotential data
        std::vector<double> p(maxLz+1);       

        for(int L=0;L<maxLz+1;++L)
        {
            double temp = utilities::BinomialFromTable(2*maxLz+2,maxLz+1);

            p[maxLz-L] = mulFactor*2.0/sqrt((double)maxLz/2.0) * 
            utilities::BinomialFromTable(2*maxLz-2*L,maxLz-L) * 
            utilities::BinomialFromTable(2*maxLz+2*L+2,maxLz+L+1)/(temp*temp);
        }
        
        return p;
    }

}   //  End namespace diagonalization 
