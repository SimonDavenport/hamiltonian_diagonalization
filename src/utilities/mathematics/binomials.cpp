////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		Contains a table of binomial numbers and a look up function. The table
//!     excludes cases where the results are either 1, N or that can be related 
//!     by N C K = N K (N-K)
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

#include "binomials.hpp" 

namespace utilities
{
    //!
    //! A function to address the binomial table (returns n choose k)
    //!
    double BinomialFromTable(
        const int n,    //!<    First binmoial argument
        const int k)    //!<    Second binomial argument
    {
        if(k==n || k==0)
        {
            return 1;
        }
        else if(n<k)
        {
            return 0;
        }
        else if(n>maxBinomial || k>maxBinomial)    //  Beyond range stored in table
        {
            return utilities::Binomial<double>(n, k);
        }
        int tmpK = k;
        if(k>n/2)
        {
            tmpK = n - tmpK;
        }
        if(tmpK==1)
        {
            return n;
        }
        else
        {
            const int tmp = n/2 - 2;
            return binomialLookUpTable[((n & 1) + tmp)*(tmp+1)+(tmpK-2)];
        }
    }
}   // End namespace utilities
