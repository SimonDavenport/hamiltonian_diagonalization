////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains template based tools to serialize small integer
//!     types into bit streams, or longer types
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

#ifndef _SERIALIZE_HPP_INCLUDED_
#define _SERIALIZE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "template_tools.hpp" 
#include <cstdint>
#include <climits>
#include <limits>
#include <sstream>
#if _DEBUG_
#include "debug.hpp"
#include <bitset>
#endif

namespace utilities
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for serializing two 8 bit values into
    //! one 16 bit value
    ////////////////////////////////////////////////////////////////////////////////
    inline uint16_t Pack2x8(
        const uint8_t key1,    //!<    Key 1 set to second 8 bits
        const uint8_t key2)    //!<    Key 2 set to first 8 bits
    {
        return ((uint16_t)key1 << 0x8) | (uint16_t)key2;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for unserializing one 16 bit value into
    //! two 8 bit values
    ////////////////////////////////////////////////////////////////////////////////
    inline void Unpack2x8(
        const uint16_t index,  //!<    Index packed with Pack2x8
        uint8_t& key1,         //!<    Key 1 set to second 8 bits
        uint8_t& key2)         //!<    Key 2 set to first 8 bits
    {
        key1 = (uint8_t)((index & 0xFF00) >> 0x8);
        key2 = (uint8_t)(index & 0x00FF);
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for serializing two 16 bit values into
    //! one 32 bit value
    ////////////////////////////////////////////////////////////////////////////////
    inline uint32_t Pack2x16(
        const uint16_t key1,    //!<    Key 1 set to second 16 bits
        const uint16_t key2)    //!<    Key 2 set to first 16 bits
    {
        return ((uint32_t)key1 << 0x10) | (uint32_t)key2;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for unserializing one 32 bit value into
    //! two 16 bit values
    ////////////////////////////////////////////////////////////////////////////////
    inline void Unpack2x16(
        const uint32_t index,   //!<    Index packed with Pack2x16
        uint16_t& key1,         //!<    Key 1 set to second 16 bits
        uint16_t& key2)         //!<    Key 2 set to first 16 bits
    {
        key1 = (uint16_t)((index & 0xFFFF0000) >> 0x10);
        key2 = (uint16_t)(index & 0x0000FFFF);
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for serializing two 32 bit values into
    //! one 64 bit value
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Pack2x32(
        const uint32_t key1,    //!<    Key 1 set to second 32 bits
        const uint32_t key2)    //!<    Key 2 set to first 32 bits
    {
        return ((uint64_t)key1 << 0x20) | (uint64_t)key2;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for unserializing one 64 bit value into
    //! two 32 bit values
    ////////////////////////////////////////////////////////////////////////////////
    inline void Unpack2x32(
        const uint64_t index,   //!<    Index packed with Pack2x32
        uint32_t& key1,         //!<    Key 1 set to second 32 bits
        uint32_t& key2)         //!<    Key 2 set to first 32 bits
    {
        key1 = (uint32_t)((index & 0xFFFFFFFF00000000) >> 0x20);
        key2 = (uint32_t)(index & 0x00000000FFFFFFFF);
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for serializing four 16 bit values into
    //! one 64 bit value
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Pack4x16(
        const uint16_t key1,    //!<    Key 1 set to fourth 16 bits
        const uint16_t key2,    //!<    Key 2 set to third 16 bits
        const uint16_t key3,    //!<    Key 3 set to second 16 bits
        const uint16_t key4)    //!<    Key 4 set to first 16 bits
    {
        const uint64_t temp1 = ((uint64_t)key1 << 0x10) | (uint64_t)key2;
        const uint64_t temp2 = ((uint64_t)key3 << 0x10) | (uint64_t)key4;
        return Pack2x32(temp1, temp2);
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for unserializing one 64 bit value into
    //! four 16 bit values
    ////////////////////////////////////////////////////////////////////////////////
    inline void Unpack4x16(
        const uint64_t index,   //!<    Index packed with Pack4x16
        uint16_t& key1,         //!<    Key 1 set to fourth 16 bits
        uint16_t& key2,         //!<    Key 2 set to third 16 bits
        uint16_t& key3,         //!<    Key 3 set to second 16 bits
        uint16_t& key4)         //!<    Key 4 set to first 16 bits
    {
        uint32_t temp1;
        uint32_t temp2;
        Unpack2x32(index, temp1, temp2);  
        Unpack2x16(temp1, key1, key2);
        Unpack2x16(temp2, key3, key4);
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for serializing eight 8 bit values into
    //! one 64 bit value
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Pack8x8(
        const uint8_t key1,     //!<    Key 1 set to eighth 8 bits
        const uint8_t key2,     //!<    Key 2 set to seventh 8 bits
        const uint8_t key3,     //!<    Key 3 set to sixth 8 bits
        const uint8_t key4,     //!<    Key 4 set to fith 8 bits
        const uint8_t key5,     //!<    Key 5 set to fourth 8 bits
        const uint8_t key6,     //!<    Key 6 set to third 8 bits
        const uint8_t key7,     //!<    Key 7 set to second 8 bits
        const uint8_t key8)     //!<    Key 8 set to first 8 bits
    {
        const uint64_t temp1 = ((uint64_t)key1 << 0x8) | (uint64_t)key2;
        const uint64_t temp2 = ((uint64_t)key3 << 0x8) | (uint64_t)key4;
        const uint64_t temp3 = ((uint64_t)key5 << 0x8) | (uint64_t)key6;
        const uint64_t temp4 = ((uint64_t)key7 << 0x8) | (uint64_t)key8;
        return Pack4x16(temp1, temp2, temp3, temp4);
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A conversion function for unserializing one 64 bit value into
    //! eight 8 bit values
    ////////////////////////////////////////////////////////////////////////////////
    inline void Unpack8x8(
        const uint64_t index,   //!<    Index packed with Pack8x8
        uint8_t& key1,          //!<    Key 1 set to eighth 8 bits
        uint8_t& key2,          //!<    Key 2 set to seventh 8 bits
        uint8_t& key3,          //!<    Key 3 set to sixth 8 bits
        uint8_t& key4,          //!<    Key 4 set to fith 8 bits
        uint8_t& key5,          //!<    Key 5 set to fourth 8 bits
        uint8_t& key6,          //!<    Key 6 set to third 8 bits
        uint8_t& key7,          //!<    Key 7 set to second 8 bits
        uint8_t& key8)          //!<    Key 8 set to first 8 bits
    {
        uint16_t temp1;
        uint16_t temp2;
        uint16_t temp3;
        uint16_t temp4;
        Unpack4x16(index, temp1, temp2, temp3, temp4);
        Unpack2x8(temp1, key1, key2);
        Unpack2x8(temp2, key3, key4);
        Unpack2x8(temp3, key5, key6);
        Unpack2x8(temp4, key7, key8);
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief g++ does not correctly implemented partial variadic template function
    //! specialization. The workaround is to use a recursive set of classes,
    //! which contain the function that we wanted to use (see below)
    ////////////////////////////////////////////////////////////////////////////////

    template <typename I,typename... Js>
    struct SerializeImpl;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to concatenate an arbitrary 
    //! number of  arguments into a single I type container. 
    //!
    //! This particular "variadic" template class is called in the case where there 
    //! are no remaining arguments to be concatenated.    
    ////////////////////////////////////////////////////////////////////////////////
    
    template<typename I>
    struct SerializeImpl<I>
    {
        static void Value(I& value, I& shift)
        {        
            return;
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to concatenate an arbitrary 
    //! number of  arguments into a single I type container. 
    //!
    //! This particular "variadic" template class is called in the case where there 
    //! remain one or more arguments to be concatenated.   
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename J, typename... Js>
    struct SerializeImpl<I,J,Js...>
    {
        static void Value(I& value, I& shift, const J& x, const Js... args)
        {
            shift -= sizeof(J)*CHAR_BIT;        //  make room to store the current value
                                                //  starting from one end of the container 
            value |= ((I)x << shift);
            utilities::SerializeImpl<I, Js...>::Value(value, shift, args...);

        }
    };
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
            
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A "variadic" template function to concatenate an arbitrary number of 
    //! arguments into a single I type container. 
    //!
    //! The variadic template argument allows for an arbitrary number of arguments, 
    //! not necessarily of the same type.   
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I,typename... Js>
    I Serialize(const Js... args)
    {
        I value = 0;
        I shift = utilities::SizeOf<Js...>()*CHAR_BIT;
        //  set the shift to the end of the container - but only go as far as required
        //  to store all of the Js type args. If we end up with "empty space" in the
        //  middle of the container then the hashing function becomes very inefficient,
        //  so instead we push any empty space to one end. 
        utilities::SerializeImpl<I,Js...>::Value(value, shift, args...);
        return value;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief g++ does not correctly implemented partial variadic template function
    //! specialization. The workaround is to use a recursive set of classes,
    //! which contain the function that we wanted to use (see below)
    ////////////////////////////////////////////////////////////////////////////////

    template <typename I, typename R, typename... Js>
    struct UnserializeArrayImpl;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to decode an arbitrary number
    //!  of arguments from a single I type container into a stringstream. 
    //!
    //! The variadic template argument allows for an arbitrary number of arguments, 
    //! not necessarily of the same type, but for writing to an array, all types
    //! should be the same to avoid allocation errors.
    ////////////////////////////////////////////////////////////////////////////////
    
    template<typename I, typename R>
    struct UnserializeArrayImpl<I, R>
    {
        static void Value(I& value, I& shift, R* returnArray)
        {        
            return;
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template function to decode an arbitrary number
    //!  of arguments from a single I type container into a stringstream. 
    //!
    //! The variadic template argument allows for an arbitrary number of arguments, 
    //! not necessarily of the same type, but for writing to an array, all types
    //! should be the same to avoid allocation errors.  
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename R, typename J, typename... Js>
    struct UnserializeArrayImpl<I, R, J, Js...>
    {
        static void Value(I& value, I& shift, R* returnArray)
        {
            shift -= sizeof(J)*CHAR_BIT;
            I mask =  ~(((I)1 << shift) - (I)1);
            I temp = value & mask;
            value ^= temp;
            if(utilities::is_same<J, R>::value)
            {
                *(returnArray) = (J)(temp >> shift);
                returnArray++;
            }
            else
            {
                std::cerr<<"\n\tERROR in UnserializeIterator(value,shift,returnArray) : returnArray type inconsistent with template arguments."<<std::endl;
                exit(EXIT_FAILURE);
            }
            utilities::UnserializeArrayImpl<I, R, Js...>::Value(value, shift, returnArray);
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
            
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A "variadic" template function to decode an arbitrary number of 
    //! arguments from a single I type container into an array.
    //!
    //! WARNING: the return array is a memory address of a pre-allocated array
    //! of type R = type J. The memory allocation must be large enough
    //! to contain Js... elements. 
    //!
    //! The variadic template argument allows for an arbitrary number of arguments, 
    //! not necessarily of the same type, but for writing to an array, all types
    //! should be the same to avoid allocation errors.     
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename R, typename... Js>
    void UnserializeArray(
        R* returnArray, //!<    array of return values EQUAL in size to the number of Js
        I& value)
    {
        I shift = utilities::SizeOf<Js...>()*CHAR_BIT;
        utilities::UnserializeArrayImpl<I, R, Js...>::Value(value, shift, returnArray);
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief g++ does not correctly implemented partial variadic template function
    //! specialization. The workaround is to use a recursive set of classes,
    //! which contain the function that we wanted to use (see below)
    ////////////////////////////////////////////////////////////////////////////////

    template <typename I, typename... Js>
    struct UnserializeImpl;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to decode an arbitrary number
    //!  of arguments from a single I type container into a stringstream. 
    //!
    //! This particular "variadic" template class is called in the case where
    //! there are no remaining arguments to be split.
    //!     
    ////////////////////////////////////////////////////////////////////////////////
    
    template <typename I>
    struct UnserializeImpl<I>
    {
        static void Value(I& value, I& shift, std::stringstream& split)
        {
            return;
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to decode an arbitrary number
    //!  of arguments from a single I type container into a stringstream. 
    //!
    //! This particular "variadic" template class is called in the case where there 
    //! remain one or more arguments to be split.
    //!
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename J, typename... Js>
    struct UnserializeImpl<I, J, Js...>
    {
        static void Value(I& value, I& shift, std::stringstream& split)
        {
            shift -= sizeof(J)*CHAR_BIT;
            I mask =  ~(((I)1 << shift) - (I)1);
            I temp = value & mask;
            value ^= temp;
            if(utilities::is_same<J, unsigned char>::value)
            {
                split<<((int)(temp >> shift))<<"\t";
            }
            else
            {
                split<<((J)(temp >> shift))<<"\t";
            }
            utilities::UnserializeImpl<I, Js...>::Value(value, shift, split);
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
            
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A "variadic" template function to decode an arbitrary number of 
    //! arguments from a single I type container into a stringstream. 
    //!
    //! The variadic template argument allows for an arbitrary number of arguments, 
    //! not necessarily of the same type.  
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename... Js>
    std::string Unserialize(
        I& value)
    {
        std::stringstream split;
        split.str(""); 
        I shift = utilities::SizeOf<Js...>()*CHAR_BIT;
        utilities::UnserializeImpl<I, Js...>::Value(value, shift, split);
        return split.str();
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief g++ does not correctly implemented partial variadic template function
    //! specialization. The workaround is to use a recursive set of classes,
    //! which contain the function that we wanted to use (see below)
    ////////////////////////////////////////////////////////////////////////////////

    template <typename I, typename... Js>
    struct SerializeFromFileImpl;
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to read in a series of parameters
    //! in a file, interpret them as a set of variadic  parameters, and finally,
    //! combine them into a single larger type (by binary concatenation)
    //!
    //! This particular "variadic" template class is called in the case where there 
    //! are no remaining arguments to be concatenated.     
    ////////////////////////////////////////////////////////////////////////////////
    
    template<typename I>
    struct SerializeFromFileImpl<I>
    {
        static void Value(std::ifstream& f_in, I& value, I& shift)
        {   
            return;
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
        
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Part of a "variadic" template class to read in a series of parameters
    //! in a file, interpret them as a set of variadic  parameters, and finally,
    //! combine them into a single larger type (by binary concatenation)
    //!
    //! This particular "variadic" template class is called in the case where there 
    //! remain one or more arguments to be concatenated.  
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename J, typename... Js>
    struct SerializeFromFileImpl<I, J, Js...>
    {
        static void Value(std::ifstream& f_in, I& value, I& shift)
        {
            shift -= sizeof(J)*CHAR_BIT;
            I temp;
            f_in>>temp;
            value |= (temp << shift);
            utilities::SerializeFromFileImpl<I, Js...>::Value(f_in, value, shift);
        }
    };
        
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
            
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A "variadic" template function to read in a series of parameters
    //! in a file, interpret them as a set of variadic  parameters, and finally,
    //! combine them into a single larger type (by binary concatenation).
    //!
    //! The variadic template argument allows for an arbitrary number of arguments, 
    //! not necessarily of the same type.  
    ////////////////////////////////////////////////////////////////////////////////

    template<typename I, typename... Js>
    I SerializeFromFile(std::ifstream& f_in)
    {
        I value = 0;
        I shift = utilities::SizeOf<Js...>()*CHAR_BIT;
        utilities::SerializeFromFileImpl<I, Js...>::Value(f_in, value, shift);
        return value;
    }
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
}   //  End namespace utilities
#endif
