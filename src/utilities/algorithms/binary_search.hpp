////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 14/03/2014
//!
//!  \file
//!		This file contains a binary search algorithm that returns the index
//!     of the element in a given array that matches a given value
//!
//!                    Copyright (C) 2014 Simon C Davenport
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

#ifndef _BINARY_SEARCH_HPP_INCLUDED_
#define _BINARY_SEARCH_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <vector>       //  For std::vector overload
#include <limits>       //  For finding the maximum value of a given integer type

#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief For a given value, determine its index with a lexicographically ordered
    //! array of values. Only the first instance is returned
    //!
    //! Binary search algorithm here taken from
    //! http://en.wikipedia.org/wiki/Binary_search_algorithm
    //!
    //! \return The first position in the array where the value occurs
    //!
    //! Return error code if unable to find the value is the maximum
    //! value of an I type given by std::numeric_limits<I>::max()
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <class T,typename I>
    I BinarySearch(
        const T& value,     //!<    Value to find
        const T* array,     //!<    Lexicographically ordered array
        I dim)              //!<    Length of array
    {
        I imax = dim-1;
        I imin = 0;
            
        #if 0   //  DEBUG

        std::cout<<"SEARCHING FOR value = "<<value<<" in array: "<<std::endl;
        
        for(int i=0;i<dim;++i)
        {
            std::cout<<array[i]<<std::endl;
        }
        
        std::cout<<std::endl;
        
        #endif

        if(value<array[imin] || value>array[imax])
        {
            return std::numeric_limits<I>::max();  //  Failed to find the value
        }

        // continue searching while [imin,imax] is not empty
        while(imax >= imin)
        {
            // calculate the midpoint for roughly equal partition
            const I imid = (imax+imin)/2;

            //std::cout<<"\n\t IMID = "<<imid<<std::endl;
            //getchar();

            const T midValue = array[imid];
            
            //std::cout<<"\n\t midValue = "<<midValue<<std::endl;

            if(value == midValue)
            {
                // key found at index imid
                //std::cout<<"RETURN: "<<imid<<std::endl;
                return imid; 
            }
            // determine which subarray to search
            else if (value > midValue)
            {
                // change min index to search upper subarray
                imin = imid + 1;
            }
            else
            {   
                //  Break out to the error code value as we've
                //  searched down to the 0th array element now
                if(imax==0) break;
              
                // change max index to search lower subarray
                imax = imid - 1;
            }
        }

        return std::numeric_limits<I>::max();  //  Failed to find the value
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Overload of the binary search function for std::vector arrays
    //!
    //! \return The first position in the array where the value occurs
    //!
    //! Return error code if unable to find the value is the maximum
    //! value of an I type given by std::numeric_limits<I>::max()
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template <class T,typename I>
    I BinarySearch(
        const T& value,                 //!<    Value to find
        const std::vector<T>& array)    //!<    Lexicographically ordered array
    {
        return utilities::BinarySearch<T,I>(value,array.data(),(I)array.size());
    }

}   //  End namespace utilities

#endif
