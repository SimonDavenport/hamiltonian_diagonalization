////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                                                  
//!	 \file
//!		This file contains some implementations of template utilities that 
//!     are included in C++ 11 's type_traits library. 
//!     
//!     Also there are some utility functions for "variadic templates".
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

#ifndef _TEMPLATE_TOOLS_HPP_INCLUDED_
#define _TEMPLATE_TOOLS_HPP_INCLUDED_

namespace utilities
{
    //!
    //! A template to check that two types are the same
    //! (as implemented by the function std::is_same in C++11)
    //!
    template<typename T, typename U>
    struct is_same 
    {
        static const bool value = false; 
    };
    
    //!
    //! A template to check that two types are the same
    //! (as implemented by the function std::is_same in C++11)
    //!
    template<typename T>
    struct is_same<T,T>  //specialization
    { 
       static const bool value = true; 
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief g++ does not correctly implemented partial variadic template function
    //! specialization. The workaround is to use a recursive set of classes,
    //! which contain the function that we wanted to use (see below)
    ////////////////////////////////////////////////////////////////////////////////
    template <typename... Ts>
    struct SizeOfImpl;
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A template class to count the total number of bytes in a variadic
    //! template parameter pack (requires c++11). Count total number of bytes - 
    //! iterator in case of a single remaining argument.
    //!
    //! (this function truncates the recursive expansion of the templates)
    ////////////////////////////////////////////////////////////////////////////////
    template <>
    struct SizeOfImpl<>
    {
        static int Value()
        {
            return 0;
        }
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A template class to count the total number of bytes in a variadic
    //! template parameter pack (requires c++11). Count total number of bytes -
    //! iterator in case of multiple parameters
    ////////////////////////////////////////////////////////////////////////////////
    template <typename T, typename... Ts>
    struct SizeOfImpl<T, Ts...>
    {
        static int Value()
        {
            return sizeof(T)+SizeOfImpl<Ts...>::Value();
        }
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A template function to count the total number of bytes in a variadic
    //! template parameter pack (requires c++11)
    ////////////////////////////////////////////////////////////////////////////////
    template <typename... Ts>
    int SizeOf()
    {
        return SizeOfImpl<Ts...>::Value();
    };
}   //  End namespace utilities
#endif
