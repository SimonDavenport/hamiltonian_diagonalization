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
    //!
    //! A conversion function for serializing two 8 bit values into
    //! one 16 bit value
    //!
    inline uint16_t Pack2x8(
        const uint8_t key1,    //!<    Key 1 set to second 8 bits
        const uint8_t key2)    //!<    Key 2 set to first 8 bits
    {
        return ((uint16_t)key1 << 0x8) | (uint16_t)key2;
    }

    //!
    //! A conversion function for unserializing one 16 bit value into
    //! two 8 bit values
    //!
    inline void Unpack2x8(
        const uint16_t index,  //!<    Index packed with Pack2x8
        uint8_t& key1,         //!<    Key 1 set to second 8 bits
        uint8_t& key2)         //!<    Key 2 set to first 8 bits
    {
        key1 = (uint8_t)((index & 0xFF00) >> 0x8);
        key2 = (uint8_t)(index & 0x00FF);
        return;
    }

    //!
    //! A conversion function for serializing two 16 bit values into
    //! one 32 bit value
    //!
    inline uint32_t Pack2x16(
        const uint16_t key1,    //!<    Key 1 set to second 16 bits
        const uint16_t key2)    //!<    Key 2 set to first 16 bits
    {
        return ((uint32_t)key1 << 0x10) | (uint32_t)key2;
    }

    //!
    //! A conversion function for unserializing one 32 bit value into
    //! two 16 bit values
    //!
    inline void Unpack2x16(
        const uint32_t index,   //!<    Index packed with Pack2x16
        uint16_t& key1,         //!<    Key 1 set to second 16 bits
        uint16_t& key2)         //!<    Key 2 set to first 16 bits
    {
        key1 = (uint16_t)((index & 0xFFFF0000) >> 0x10);
        key2 = (uint16_t)(index & 0x0000FFFF);
        return;
    }

    //!
    //! A conversion function for serializing two 32 bit values into
    //! one 64 bit value
    //!
    inline uint64_t Pack2x32(
        const uint32_t key1,    //!<    Key 1 set to second 32 bits
        const uint32_t key2)    //!<    Key 2 set to first 32 bits
    {
        return ((uint64_t)key1 << 0x20) | (uint64_t)key2;
    }

    //!
    //! A conversion function for unserializing one 64 bit value into
    //! two 32 bit values
    //!
    inline void Unpack2x32(
        const uint64_t index,   //!<    Index packed with Pack2x32
        uint32_t& key1,         //!<    Key 1 set to second 32 bits
        uint32_t& key2)         //!<    Key 2 set to first 32 bits
    {
        key1 = (uint32_t)((index & 0xFFFFFFFF00000000) >> 0x20);
        key2 = (uint32_t)(index & 0x00000000FFFFFFFF);
        return;
    }

    //!
    //! A conversion function for serializing four 16 bit values into
    //! one 64 bit value
    //!
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

    //!
    //! A conversion function for unserializing one 64 bit value into
    //! four 16 bit values
    //!
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

    //!
    //! A conversion function for serializing eight 8 bit values into
    //! one 64 bit value
    //!
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

    //!
    //!  A conversion function for unserializing one 64 bit value into
    //! eight 8 bit values
    //!
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
}   //  End namespace utilities
#endif
