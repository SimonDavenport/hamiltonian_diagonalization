////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This file contains a class to implement a multi-key hash table and a
//!     class to implement a multi-key multi-value hash table. 
//!
//!     The class is a wrapper for an underlying single-key hash table 
//!     implementations defined by STL, or by "sparsehash", where available.
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

#ifndef _MULTI_KEY_HASH_HPP_INCLUDED_
#define _MULTI_KEY_HASH_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../wrappers/mpi_wrapper.hpp"
#include "../wrappers/io_wrapper.hpp"
#include "../wrappers/murmur_hash_wrapper.hpp"
#include "../general/serialize.hpp"
#include "../../utilities/algorithms/quick_sort.hpp"
#if _SPEED_OPTIMIZED_MAP_
#include <sparsehash/dense_hash_map>
#elif _MEMORY_OPTIMIZED_MAP_
#include <sparsehash/sparse_hash_map>
#endif
#include <unordered_map>
#include <vector>
#include <fstream>
#include <set>
#if _DEBUG_
#include "../general/debug.hpp"
#include <bitset>
#endif

namespace utilities
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Function to generate a 64-bit hash key from 2 input arguments
    //!
    //! \return combined 64-bit key
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Key(
        const uint32_t key1,    //!<    Key 1
        const uint32_t key2)    //!<    Key 2
    {
        return utilities::Pack2x32(key1, key2);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Function to generate a 64-bit hash key from 3 input arguments
    //!
    //! \return combined 64-bit key
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Key(
        const uint16_t key1,    //!<    Key 1
        const uint16_t key2,    //!<    Key 2
        const uint16_t key3)    //!<    Key 3
    {
        return utilities::Pack4x16(key1, key2, key3, 0);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Function to generate a 64-bit hash key from 4 input arguments
    //!
    //! \return combined 64-bit key
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Key(
        const uint16_t key1,    //!<    Key 1
        const uint16_t key2,    //!<    Key 2
        const uint16_t key3,    //!<    Key 3
        const uint16_t key4)    //!<    Key 4
    {
        return utilities::Pack4x16(key1, key2, key3, key4);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Function to generate a 64-bit hash key from 6 input arguments
    //!
    //! \return combined 64-bit key
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Key(
        const uint8_t key1,    //!<    Key 1
        const uint8_t key2,    //!<    Key 2
        const uint8_t key3,    //!<    Key 3
        const uint8_t key4,    //!<    Key 4
        const uint8_t key5,    //!<    Key 5
        const uint8_t key6)    //!<    Key 6
    {
        return utilities::Pack8x8(key1, key2, key3, key4, key5, key6, 0, 0);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Function to generate a 64-bit hash key from 8 input arguments
    //!
    //! \return combined 64-bit key
    ////////////////////////////////////////////////////////////////////////////////
    inline uint64_t Key(
        const uint8_t key1,     //!<    Key 1
        const uint8_t key2,     //!<    Key 2
        const uint8_t key3,     //!<    Key 3
        const uint8_t key4,     //!<    Key 4
        const uint8_t key5,     //!<    Key 5
        const uint8_t key6,     //!<    Key 6
        const uint8_t key7,     //!<    Key 7
        const uint8_t key8)     //!<    Key 8
    {
        return utilities::Pack8x8(key1, key2, key3, key4, key5, key6, key7, key8);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A class template to contain implementations of utility 
    //! functions used in either the MultiHashMap or MultiHashMultiMap 
    //! classes below
    //!
    //! The basic multi hash map is built by a binary concatenation of  
    //! multiple keys, which are then treated as a single value for an 
    //! underlying hash  table implementation. 
    //!
    //! The template parameter T defines the type of data to be stored   
    //! in the hash map
    ////////////////////////////////////////////////////////////////////////////////   
    template<typename T>
    class MultiHashBase
    {
        protected:
        ////////////////////////////////////////////////////////////////////////////////
        #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to return the default deleted map key
        //!
        //! Currently set to be the highest possible value 
        //! of the template integer type
        ////////////////////////////////////////////////////////////////////////////////
        uint64_t DeletedMapKey() const
        {
            return std::numeric_limits<uint64_t>::max();
        }
        #endif
        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////
        #if _SPEED_OPTIMIZED_MAP_ 
        
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief A function to return the default empty map key
        //!
        //! Currently set to be the highest possible value 
        //! of the template integer type minus 1
        ////////////////////////////////////////////////////////////////////////////////
        uint64_t EmptyMapKey() const
        {
            return std::numeric_limits<uint64_t>::max()-1;
        }
        #endif
        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Mpi synchronize the contents of the map. This template function
        //! takes any map type as an argument, but every map type has available all
        //! the methods called for this operation.
        ////////////////////////////////////////////////////////////////////////////////
        template <class M>
        void MpiSynchronizeBase(
            M& map,                 //!<    Map to be synchronized
            const int nodeId,       //!<    Node to synchronize with
            const MpiWrapper& mpi)  //!<    Instance of the MPI wrapper class
        {
            //  This implementation is very inefficient: we simply read out the 
            //  contents of the map into a pair of key and map arrays, send
            //  these arrays over the network, then reconstruct the hash tables
            //  on the other nodes
            //  First broadcast the number of elements in the map
            uint64_t size;
            if(nodeId == mpi.m_id)    //  FOR THE MASTER NODE
            {
                size = map.size();
            }
            MPI_Bcast(&size, 1, mpi.GetType<uint64_t>(), nodeId, mpi.m_comm);
            //  Allocate buffers to store the map on each node
            uint64_t* keyBuffer   = new uint64_t[size];
            T* valueBuffer = new T[size];
            if(nodeId == mpi.m_id)    //  FOR THE MASTER NODE
            {
                //  Populate the buffers on the master node
                uint64_t* p_keyBuffer = keyBuffer;
                T* p_valueBuffer = valueBuffer;
                for(auto it = map.begin(); it != map.end(); ++it, ++p_keyBuffer, ++p_valueBuffer)
                {
                    *(p_keyBuffer) = it->first;
                    *(p_valueBuffer) = it->second;
                }
            }
            MPI_Bcast(keyBuffer, size, mpi.GetType<uint64_t>(), nodeId, mpi.m_comm); 
            MPI_Bcast(valueBuffer, size,mpi.GetType<T>(), nodeId, mpi.m_comm);
            //  Set the maps on all other nodes from the buffered values
            //  (the map values on other nodes are cleared first)
            if(nodeId != mpi.m_id)    //  FOR ALL OTHER NDOES
            {
                if(!map.empty())
                {
                    map.clear();
                }
                uint64_t* p_keyBuffer = keyBuffer;
                T* p_valueBuffer = valueBuffer;
                //  Generate a pair buffer so that all map values can 
                //  be inserted in one operation
                std::vector<std::pair<uint64_t, T> > pairList;
                pairList.reserve(size);
                for(uint64_t i=0; i<size; ++i, ++p_keyBuffer, ++p_valueBuffer)
                {
                    pairList.push_back(std::pair<uint64_t, T>(*(p_keyBuffer), *(p_valueBuffer)));
                }
                //  Insert the whole range into the map in one operation
                map.insert(pairList.begin(), pairList.end());
            }
            delete[] keyBuffer;
            delete[] valueBuffer;
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Mpi gather the contents of the map onto the master node. 
        //! This template function takes any map type as an argument, but every 
        //! map type has available all the methods called for this operation.
        ////////////////////////////////////////////////////////////////////////////////
        template <class M>
        void MpiGatherBase(
            M& map,                 //!<    Map to be gathered
            const int nodeId,       //!<    Node to gather on
            const MpiWrapper& mpi)  //!<    Instance of the MPI wrapper class
        {
            //  First count the number of values that are in the map on each node
            int size;
            int totalSize;
            size = map.size();
            MPI_Reduce(&size, &totalSize, 1, mpi.GetType<int>(), MPI_SUM, nodeId, mpi.m_comm);
            //  Convert all other node maps to a list of buffers for the keys and values
            uint64_t* keyBuffer = new uint64_t[size];
            T* valueBuffer      = new T[size];
            uint64_t* p_keyBuffer = keyBuffer;
            T* p_valueBuffer      = valueBuffer;
            for(auto it = map.begin(); it != map.end(); ++it, ++p_keyBuffer, ++p_valueBuffer)
            {
                *(p_keyBuffer)   = it->first;
                *(p_valueBuffer) = it->second;
            }
            //  Clear the maps on all nodes to free up available space
            if(!map.empty())
            {
                map.clear();
            }
            //  Allocate memory to store the accumulated buffers on the master node
            uint64_t* keyRecvBuffer = 0;
            T* mapRecvBuffer = 0;
            if(nodeId == mpi.m_id)    //  For the gather node
            {
                keyRecvBuffer = new uint64_t[totalSize];
                mapRecvBuffer = new T[totalSize];
            }
            //  Use an MPI gather routine to move around the buffered data
            MPI_Status status;
            mpi.Gather<uint64_t>(keyBuffer, size, keyRecvBuffer, totalSize, nodeId, mpi.m_comm, status);
            mpi.Gather<T>(valueBuffer, size, mapRecvBuffer, totalSize, nodeId, mpi.m_comm, status);
            delete[] keyBuffer;
            delete[] valueBuffer;
            //  Insert the received buffer data into a map on the master node
            if(nodeId == mpi.m_id)    //  For the gather node
            {
                //  First convert the buffered map and value arrays into 
                //  a pair list
                std::vector<std::pair<uint64_t, T> > pairList;
                pairList.reserve(totalSize);
                uint64_t* p_keyBuffer = keyRecvBuffer;
                T* p_valueBuffer = mapRecvBuffer;
                for(uint64_t i=0; i<totalSize; ++i, ++p_keyBuffer, ++p_valueBuffer)
                {
                    pairList.push_back(std::pair<uint64_t, T>(*(p_keyBuffer), *(p_valueBuffer)));
                }
                //  Remove duplicates
                QuickSort<std::pair<uint64_t, T>, int, _ASCENDING_ORDER_>(pairList);
                std::vector<std::pair<uint64_t, T> > uniquePairList;
                uniquePairList.reserve(totalSize);
                std::pair<uint64_t, T> prev_value = pairList[0];
                for(typename std::vector<std::pair<uint64_t, T> >::const_iterator it = pairList.begin(); 
                    it < pairList.end(); ++it)
                {
                    uniquePairList.push_back(*it);
                    while(*it == prev_value)
                    {
                        ++it;
                    }
                    prev_value = *it;
                }
                //  Insert the whole range into the map in one operation
                map.insert(uniquePairList.begin(), uniquePairList.end());
                delete[] keyRecvBuffer;
                delete[] mapRecvBuffer;  
            }
            return;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Mpi scatter the contents of the map from the master node to
        //! distribute over a number of nodes according to the range of the map
        //! key 
        //!
        //! This template function takes any map type as an argument, but every 
        //! map type has available all the methods called for this operation.
        ////////////////////////////////////////////////////////////////////////////////
        template <class M>
        void MpiScatterBase(
            M& map,                 //!<    Map to be gathered
            const uint64_t maxKey,  //!<    Maximum map key stored on the current node
            const int nodeId,       //!<    Node to scatter from
            MpiWrapper& mpi)        //!<    Instance of the MPI wrapper class
        {
            //  Communicate all of the maxKey values to every node
            uint64_t maxKeyList[mpi.m_nbrProcs];
            maxKeyList[mpi.m_id] = maxKey;
            for(unsigned int p=0; p<mpi.m_nbrProcs; ++p)
            {
                mpi.Sync(&maxKeyList[p], 1, p);
            }
            //  Check that the maxKeyList is in ascending order
            for(unsigned int p=1; p<mpi.m_nbrProcs; ++p)
            {
                if(maxKeyList[p]<maxKeyList[p-1])
                {
                    if(nodeId == mpi.m_id)    //  For the scatter node
                    {
                        std::cerr<<"\n\tERROR in MpiScatter: maxKey declaration not in asceding order!"<<std::endl;
                    }
                    exit(EXIT_FAILURE);
                }
            }
            //  Prepare a list of the sizes of the hash map on each node
            unsigned int sizeList[mpi.m_nbrProcs];
            for(unsigned int p=0; p<mpi.m_nbrProcs; ++p)
            {
                sizeList[p] = 0;
            }
            uint64_t* keyBuffer = 0;
            T* valueBuffer = 0;
            unsigned int totalSize;
            if(nodeId == mpi.m_id)    //  For the scatter node
            {
                //  First count the number of values that are in the map                
                totalSize = map.size();
                //  Convert all other node maps to a list of buffers for the keys and values
                keyBuffer   = new uint64_t[totalSize];
                valueBuffer = new T[totalSize];
                uint64_t* p_keyBuffer   = keyBuffer;
                T* p_valueBuffer = valueBuffer;
                for(auto it = map.begin(); it != map.end(); ++it, ++p_keyBuffer, ++p_valueBuffer)
                {
                    *p_keyBuffer   = it->first;
                    *p_valueBuffer = it->second;
                }
                //  Clear the map to free up available space
                if(!map.empty())
                {
                    map.clear();
                    #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
                    map.resize(0);
                    #endif
                }
                //  Count how many map key values are going to be stored on each
                //  node once the map is distributed
                p_keyBuffer = keyBuffer;

                for(unsigned int i=0; i<totalSize; ++i, ++p_keyBuffer)
                {
                    for(unsigned int p=0; p<mpi.m_nbrProcs; ++p)
                    {
                        if(*p_keyBuffer < maxKeyList[p])
                        {
                            ++sizeList[p];
                            break;
                        }
                    }
                }
                unsigned int cumulativeSizeList[mpi.m_nbrProcs];
                for(unsigned int p=0; p<mpi.m_nbrProcs; ++p)
                {
                    cumulativeSizeList[p] = 0;
                    for(unsigned int q=0; q<p; ++q)
                    {
                        cumulativeSizeList[p] += sizeList[q];
                    }
                }
                //  Check that the total of the size list equals the size of
                //  the hash map
                if(cumulativeSizeList[mpi.m_nbrProcs-1]+sizeList[mpi.m_nbrProcs-1] != totalSize)
                {
                    std::cerr<<"\n\tERROR in MultiHashMap::MpiScatter: list of maximum keys leads to inconsistent data distribution."<<std::endl;
                    mpi.m_exitFlag = true;
                }
                //  Reorder the buffers so that the data is in the correct order
                //  to be contiguously distributed to the other nodes
                uint64_t* reorderKeyBuffer = new uint64_t[totalSize];
                T* reorderValueBuffer = new T[totalSize];
                p_keyBuffer = keyBuffer;
                p_valueBuffer = valueBuffer;
                for(unsigned int p=0; p<mpi.m_nbrProcs; ++p)
                {
                    sizeList[p] = 0;
                }
                for(unsigned int i=0; i<totalSize; ++i, ++p_keyBuffer, ++p_valueBuffer)
                {
                    for(unsigned int p=0; p<mpi.m_nbrProcs; ++p)
                    {
                        if(*p_keyBuffer < maxKeyList[p])
                        {
                            reorderKeyBuffer[cumulativeSizeList[p]+sizeList[p]] = *p_keyBuffer;
                            reorderValueBuffer[cumulativeSizeList[p]+sizeList[p]] = *p_valueBuffer;
                            ++sizeList[p];
                            break;
                        }
                    }
                }
                memcpy(keyBuffer, reorderKeyBuffer, totalSize*sizeof(uint64_t));
                memcpy(valueBuffer, reorderValueBuffer, totalSize*sizeof(T));
                delete[] reorderKeyBuffer;
                delete[] reorderValueBuffer;
            }
            mpi.ExitFlagTest();
            //  Next, inform each of the other nodes about the number of values
            //  they should receive
            mpi.Sync(sizeList,mpi.m_nbrProcs, nodeId);
            mpi.Sync(&totalSize, 1, 0);
            uint64_t* keyRecvBuffer = 0;
            T* valueRecvBuffer = 0;
            keyRecvBuffer = new uint64_t[sizeList[mpi.m_id]];
            valueRecvBuffer = new T[sizeList[mpi.m_id]];
            //  Scatter the key and values buffers over all nodes
            MPI_Status status;
            mpi.Scatter(keyBuffer, totalSize, keyRecvBuffer, sizeList[mpi.m_id], nodeId, mpi.m_comm, status);
            mpi.Scatter(valueBuffer, totalSize, valueRecvBuffer, sizeList[mpi.m_id], nodeId, mpi.m_comm, status);
            if(0 != keyBuffer)   delete[] keyBuffer;
            if(0 != valueBuffer) delete[] valueBuffer;
            //  Finally, insert the key and map back into a hash table object
            //  and then delete the buffers
            std::vector<std::pair<uint64_t, T> > pairList;
            pairList.reserve(sizeList[mpi.m_id]);
            uint64_t* p_keyBuffer = keyRecvBuffer;
            T* p_valueBuffer = valueRecvBuffer;
            for(uint64_t i=0; i<sizeList[mpi.m_id]; ++i, ++p_keyBuffer, ++p_valueBuffer)
            {
                pairList.push_back(std::pair<uint64_t, T>(*(p_keyBuffer), *(p_valueBuffer)));
            }
            //  Insert the whole range into the map in one operation
            map.insert(pairList.begin(), pairList.end());
            //  Finally remove memory allocated to store the receive buffers
            delete[] keyRecvBuffer;
            delete[] valueRecvBuffer;
        }

        //!
        //! Write the contents of the map to a file using the set number
        //! of key labels
        //!
        template <class M>
        void ToFileBase(
            const M& map,                   //!<    Map containing data to write to the file
            std::ofstream& stream,          //!<    File stream
            const io::fileFormat_t format,  //!<    Format of file
            const unsigned int nbrLabels,   //!<    Number of key labels
            MpiWrapper& mpi)                //!<    Instance of the MPI wrapper class  
        {
            //  Create a copy of the map for the MPI gather operation to modify
            M temp_map = map;
            this->MpiGatherBase(temp_map, 0, mpi);
            if(0 == mpi.m_id)    //  FOR THE MASTER NODE
            {
                if(!mpi.m_exitFlag)
                {
                    if(io::_BINARY_ == format)
                    {   
                        unsigned int size = temp_map.size();
                        std::vector<uint64_t> keys(size);
                        std::vector<T> values(size);
                        std::vector<uint64_t>::iterator keys_it = keys.begin();
                        typename std::vector<T>::iterator values_it = values.begin();
                        for(auto it = temp_map.begin(); it != temp_map.end(); ++it, ++keys_it, ++values_it)
                        {
                            *keys_it = it->first;
                            *values_it = it->second;
                        }
                        stream.write((char*)&size, sizeof(unsigned int));
                        stream.write((char*)keys.data(), size*sizeof(uint64_t));
                        stream.write((char*)values.data(), size*sizeof(T));
                    }
                    else if(io::_TEXT_ == format)
                    {
                        stream << temp_map.size() << "\n";
                        stream << nbrLabels << "\n";
                        //  Write map key-value pairs to the file (including multiple  
                        //  values if found)
                        for(auto it = temp_map.begin(); it != temp_map.end(); ++it)
                        {
                            std::string streamData = utilities::ToStream(it->second);
                            switch(nbrLabels)
                            {
                                case 1:
                                {
                                    stream << it->first << "\t" << streamData <<"\n";
                                    break;
                                }
                                case 2:
                                {    
                                    uint32_t key1, key2;
                                    utilities::Unpack2x32(it->first, key1, key2);
                                    stream << key1 << "\t" << key2 << "\t" << streamData <<"\n";
                                    break;
                                }
                                case 3:
                                {    
                                    uint16_t key1, key2, key3, key4;
                                    utilities::Unpack4x16(it->first, key1, key2, key3, key4);
                                    stream << key1 << "\t" << key2 << "\t" << key3 << "\t" << streamData <<"\n";
                                    break;
                                }
                                case 4:
                                {
                                    uint16_t key1, key2, key3, key4;
                                    utilities::Unpack4x16(it->first, key1, key2, key3, key4);           
                                    stream << key1 << "\t" << key2 << "\t" << key3 << "\t" << key4 << "\t";
                                    stream << streamData <<"\n";
                                    break;
                                }
                                case 8:
                                {
                                    uint8_t key1, key2, key3, key4, key5, key6, key7, key8;
                                    utilities::Unpack8x8(it->first, key1, key2, key3, key4, key5, key6, key7, key8);
                                    stream << key1 << "\t" << key2 << "\t" << key3 << "\t" << key4 << "\t" << key5<< "\t" << key6;
                                    stream << "\t" << key7 << "\t" << key8 << "\t" << streamData <<"\n";
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            mpi.ExitFlagTest();
            return;
        }

        //!
        //! Write the contents of the map to a file using the set number
        //! of key labels
        //!
        template <class M>
        void FromFileBase(
            M& map,                         //!<    Map containing data to write to the file
            std::ifstream& stream,          //!<    File stream
            const io::fileFormat_t format,  //!<    Format of file
            MpiWrapper& mpi)                //!<    Instance of the MPI wrapper class  
        {
            if(0 == mpi.m_id)    //  FOR THE MASTER NODE
            {
                if(!mpi.m_exitFlag)
                {
                    std::vector<std::pair<uint64_t, T> > pairList;
                    if(io::_BINARY_ == format)
                    {
                        unsigned int size = 0;
                        stream.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
                        pairList.resize(size);
                        std::vector<uint64_t> keys(size);
                        std::vector<T> values(size);
                        stream.read(reinterpret_cast<char*>(keys.data()), size*sizeof(uint64_t));
                        stream.read(reinterpret_cast<char*>(values.data()), size*sizeof(T));
                        std::vector<uint64_t>::iterator keys_it = keys.begin();
                        typename std::vector<T>::iterator values_it = values.begin();
                        for(auto it = pairList.begin(); it < pairList.end(); ++it, ++keys_it, ++values_it)
                        {
                            it->first = *keys_it;
                            it->second = *values_it;
                        }
                    }
                    else if(io::_TEXT_ == format)
                    {
                        unsigned int size;
                        stream >> size;
                        uint64_t nbrLabels; 
                        stream >> nbrLabels;
                        pairList.reserve(size);
                        for(uint64_t i=0; i<size; ++i)
                        {
                            uint64_t key;
                            switch(nbrLabels)
                            {
                                case 1:
                                {
                                    stream >> key;
                                    break;
                                }
                                case 2:
                                {
                                    uint32_t key1, key2;
                                    stream >> key1;
                                    stream >> key2;
                                    key = utilities::Pack2x32(key1, key2);
                                    break;
                                }
                                case 3:
                                {
                                    uint16_t key1, key2, key3;
                                    stream >> key1;
                                    stream >> key2;
                                    stream >> key3;
                                    key = utilities::Pack4x16(key1, key2, key3, 0);
                                    break;
                                }
                                case 4:
                                {
                                    uint16_t key1, key2, key3, key4;
                                    stream >> key1;
                                    stream >> key2;
                                    stream >> key3;
                                    stream >> key4;
                                    key = utilities::Pack4x16(key1, key2, key3, key4);
                                    break;
                                }
                                case 8:
                                {
                                    uint8_t key1, key2, key3, key4, key5, key6, key7, key8;
                                    stream >> key1;
                                    stream >> key2;
                                    stream >> key3;
                                    stream >> key4;
                                    stream >> key5;
                                    stream >> key6;
                                    stream >> key7;
                                    stream >> key8;
                                    key = utilities::Pack8x8(key1, key2, key3, key4, key5, key6, key7, key8);
                                    break;
                                }
                            }
                            T val;
                            utilities::FromStream(stream, val);
                            pairList.push_back(std::pair<uint64_t, T>(key, val));
                        }
                    }
                    map.insert(pairList.begin(), pairList.end());
                }
            }
            mpi.ExitFlagTest();
            this->MpiSynchronizeBase(map, 0, mpi);
            return;
        }
        public:
        
        //!
        //! MultiHashBase constructor
        //!
        MultiHashBase(){};
        
        //!
        //! MultiHashBase destructor
        //!
        ~MultiHashBase(){};
    }; 

    //!
    //! A class template to implement a multi-key hash table. 
    //! 
    template<typename T>
    class MultiHashMap : MultiHashBase<T>
    {
        private:
        T m_emptyReturn;                        //!<    Value to return in case 
                                                //!     of no map entry found
        public:
        ////////////////////////////////////////////////////////////////////////////////
        #if _SPEED_OPTIMIZED_MAP_
        google::dense_hash_map<uint64_t, T, MurmurHasher64Wrapper<uint64_t>> m_map;     
                                                //!<    An extremely fast hash map
                                                //!     implementation          
        #elif _MEMORY_OPTIMIZED_MAP_
        google::sparse_hash_map<uint64_t, T, MurmurHasher64Wrapper<uint64_t>> m_map;     
                                                //!<    A memory efficient hash table
                                                //!     container to represent the matrix
        #else
        std::unordered_map<uint64_t, T> m_map;  //!<     Default to using the STL
                                                //!      unordered map
        #endif
        ////////////////////////////////////////////////////////////////////////////////
          
        //!
        //! Default constructor sets empty map key and deleted map key
        //!
        MultiHashMap()
        :   MultiHashBase<T>(),
            m_emptyReturn(0)
        {
            #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
            m_map.set_deleted_key(this->DeletedMapKey());
            #endif
            #if _SPEED_OPTIMIZED_MAP_
            m_map.set_empty_key(this->EmptyMapKey());
            #endif
            return;
        }

        //!
        //! Destructor clears the hash map implementation
        //!
        ~MultiHashMap()
        {   
            this->Clear();
        }
          
        //!
        //! Define a copy assignment operator
        //!        
        MultiHashMap<T>& operator = (const MultiHashMap<T>& other)
        {
            //  First erase the currently stored map
            this->Clear();
            //  Then make a copy assignment to another map
            m_map = other.m_map;
            return *this;
        }
           
        //!
        //! The clear function clears the current hash map.
        //!
        void Clear()
        {
            if(!m_map.empty())
            {
                m_map.clear();
            }
            #if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
            m_map.resize(0);
            #endif
        }
        
        //!
        //! Set a new value for the error code
        //!
        void SetEmptyReturn(
            const T emptyReturn)        //!<    Value to return in case of no 
                                        //!     map entry found
        {
            m_emptyReturn = emptyReturn;
        }

        ////////////////////////////////////////////////////////////////////////////////
        //! \brief Function to append data to the map using a single index
        //!
        //! \return The address of the matrix element, which can then be set with 
        //! an = operation e.g. Insert(i) = value
        ////////////////////////////////////////////////////////////////////////////////
        T& Insert(
            const uint64_t key1)    //!<    Key 1
        {
            return m_map[key1];
        }

        //!
        //! Function to retrun address of map values using a single index
        //! 
        const T Value(
            const uint64_t index)   //!<    Map index
            const
        {
            auto it = m_map.find(index);
            if(m_map.end() == it)
            {
                //  Index not found in map, so return m_emptyReturn
                return m_emptyReturn;
            }
            else
            {
                return it->second;
            }
        }

        //!
        //! Function to return the number of set elements
        //!
        uint64_t Size() const
        {
            return m_map.size();
        }

        //!
        //! Store the map in a file
        //! 
        void ToFile(
            std::ofstream& stream,          //!<    File stream
            const io::fileFormat_t format,  //!<    Format of file
            const unsigned int nbrLabels,   //!<    Number of map key labels
            MpiWrapper& mpi)                //!<    Instance of the MPI wrapper class
        {
            this->ToFileBase(m_map, stream, format, nbrLabels, mpi);
        }

        //!
        //! Retrieve map data from a file and use it to construct a new map
        //! 
        void FromFile(
            std::ifstream& stream,          //!<    File stream
            const io::fileFormat_t format,  //!<    Format of file
            MpiWrapper& mpi)                //!<    Instance of the MPI wrapper class
        {
            this->Clear();
            this->FromFileBase(m_map, stream, format, mpi);
        }   

        //!
        //! Mpi synchronize the contents of the map
        //!
        void MpiSynchronize(
            const int nodeId,           //!<    Node to synchronize with
            const MpiWrapper& mpi)      //!<    Instance of the MPI wrapper class
        {
            this->MpiSynchronizeBase(m_map, nodeId, mpi);
        }

        //!
        //! Mpi gather the contents of parallel distributed maps into a single
        //! map on a specified node. 
        //!
        void MpiGather(
            const int nodeId,           //!<    Node to gather on
            const MpiWrapper& mpi)      //!<    Instance of the MPI wrapper class
        {
            this->MpiGatherBase(m_map, nodeId, mpi);
        }
        
        //!
        //! Mpi scatter the contents of a single map to several nodes, based on 
        //! a list of maximum map keys allowed on each node (in ascending order)
        //!   
        void MpiScatter(
            const uint64_t maxKey,      //!<    Maximum map key stored on the current node
            const int nodeId,           //!<    Node to scatter from
            const MpiWrapper& mpi)      //!<    Instance of the MPI wrapper class
        {
            this->MpiScatterBase(m_map, maxKey, nodeId, mpi);
        }
    };  //  End MultiHashMap

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //!
    //! A class template to implement a multi-key multi-hash table. 
    //!
    template<typename T>
    class MultiHashMultiMap : MultiHashBase<T>
    {
        public:
        std::unordered_multimap<uint64_t, T, MurmurHasher64Wrapper<uint64_t>> m_map;          
            //!<    Use the STL unordered_multimap
   
        //!
        //! Default constructor declared empty
        //!
        MultiHashMultiMap()
        :   MultiHashBase<T>()
        {}
           
        //!
        //! Destructor clears the hash map implementation
        //!
        ~MultiHashMultiMap()
        {
            this->Clear();
        }
            
        //!
        //! Define a copy assignment operator
        //!
        MultiHashMultiMap<T>& operator = (const MultiHashMultiMap<T>& other)
        {
            //  First erase the currently stored map
            this->Clear();
            //  Then make a copy assignment to another map
            m_map = other.m_map;
            return *this;
        }
        
        //!
        //! The clear function clears the current hash map
        //!
        void Clear()
        {
            if(!m_map.empty())
            {
                m_map.clear();
            }
        }

        //!
        //! Function to append data to the map using a single index.
        //! Use as Insert(value,i). More than one value can be subsequently
        //! inserted for the same key i. 
        //!
        void Insert(
            const T& element,           //!<  Value to be inserted
            const uint64_t key)         //!<  Map key
        {
            m_map.insert(std::pair<uint64_t, T>(key, element));
        }

        //!
        //! \brief Function to count the number of entries for a given map key
        //!
        //!
        uint32_t Count(
            const uint64_t key)     //!<  Map key
            const
        {
            return m_map.count(key);
        }
        
        //!
        //! Function to return the highese count of values stored for any 
        //! given key (can be used to allocate a bugger to store Value returns)
        //!
        uint32_t GetMaxCount()
            const
        {
            uint32_t max = 0;
            for(auto it = m_map.begin(); it!=m_map.end(); ++it)
            {
                max = std::max(max, (uint32_t)m_map.count(it->first)); 
            }
            return max;
        }

        //!
        //! Function to access values in the map. In the multi-map case this 
        //! function returns a vector of all values in the map corresponding to the 
        //! specified key.
        //!
        void Value(
            T* value,               //!<    A buffer that is large enough to store the 
                                    //!     maximum number of values for a given key
                                    //!     (use GetMaxCount() function to allocate)
            uint32_t& nbr,          //!<    Address to set actual size of return buffer
                                    //!     used in the current call           
            const uint64_t key)     //!<    Map key
            const
        {
            //  Find all values corresponding to the specified key
            auto range = m_map.equal_range(key);
            if(m_map.end() == range.first)
            {
                //  Index not found in map
                return;
            }
            T* p_value = value;
            nbr = 0;
            for (auto it=range.first; it!=range.second; ++it, ++p_value, ++nbr)
            {
                *(p_value) = it->second;
            }
            return;
        }

        //!
        //! Store the map in a file
        //! 
        void ToFile(
            std::ofstream& stream,          //!<    File stream
            const io::fileFormat_t format,  //!<    Format of file
            const unsigned int nbrLabels,   //!<    Number of map key labels
            MpiWrapper& mpi)                //!<    Instance of the MPI wrapper class
        {
            this->ToFileBase(m_map, stream, format, nbrLabels, mpi);
        }

        //!
        //! Retrieve map data from a file and use it to construct a new map
        //!
        void FromFile(
            std::ifstream& stream,          //!<    File stream
            const io::fileFormat_t format,  //!<    Format of file
            MpiWrapper& mpi)                //!<    Instance of the MPI wrapper class
        {
            this->Clear();
            this->FromFileBase(m_map, stream, format, mpi);
        }

        //!
        //! Mpi synchronize the contents of the map
        //!
        void MpiSynchronize(
            const int nodeId,           //!<    Node to synchronize with
            const MpiWrapper& mpi)      //!<    Instance of the MPI wrapper class
        {
            this->MpiSynchronizeBase(m_map, nodeId, mpi);
        }

        //!
        //! Mpi gather the contents of parallel distributed maps into a single
        //! map on a specified node. 
        //!
        void MpiGather(
            const int nodeId,           //!<    Node to gather on
            const MpiWrapper& mpi)      //!<    Instance of the MPI wrapper class
        {
            this->MpiGatherBase(m_map, nodeId, mpi);
        }
    };  //  End class MultiHashMultiMap
}  //  End namespace utilities
#endif
