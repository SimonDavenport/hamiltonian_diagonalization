////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a wrapper for the MurmurHash2 and MurmurHash3 hash 
//!     functions. The source code is available here:
//!     https://sites.google.com/site/murmurhash/
//!
//!     svn checkout http://smhasher.googlecode.com/svn/trunk/ smhasher-read-only
//!
//!     I have added some additional modifications to interface with
//!     SparseHash 
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

#ifndef _MURMUR_HASH_WRAPPER_HPP_INCLUDED_
#define _MURMUR_HASH_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../mathematics/murmur_hash2.hpp"
#if _DEBUG_
#include "../general/debug.hpp"
#endif

namespace utilities
{
    //!
    //! Template wrapper class for interfacing with SparseHash
    //!
    template<typename T>
    struct MurmurHasher32Wrapper
    {
        uint32_t operator()(const T& t) const 
        {
            return MurmurHash2A(&t, sizeof(t), 0);
        }
    };

    //!
    //! Template wrapper class for interfacing with SparseHash
    //!
    template<typename T>
    struct MurmurHasher64Wrapper
    {
        uint64_t operator()(const T& t) const 
        {
            return MurmurHash64A(&t, sizeof(t), 0);
        }
    };
    
    //!
    //! Template wrapper function for general use of seeded hash functions
    //!
    template<typename T>
    uint64_t MurmurHasher64(const T& t, const unsigned int seed) 
    {
        return MurmurHash64A(&t, sizeof(t), seed);
    }
}   //  End namespace utilities
#endif
