////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store pre-calculated Clebsch-Gordan 
//!     coefficients. Coefficients are calculated using the algorithm described 
//!     here: arXiv:1009.0437v1. See:
//!     http://homepages.physik.uni-muenchen.de/~vondelft/Papers/ClebschGordan/
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

#ifndef _CLEBSCH_GORDAN_COEFFICIENTS_HPP_INCLUDED_
#define _CLEBSCH_GORDAN_COEFFICIENTS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../wrappers/lapack_wrapper.hpp"
#include <cstdint>
#include <cmath>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <vector>
#include <map>
#include <functional>       
#if _DEBUG_
#include "../general/debug.hpp"
#endif

//////////////////////////////////////////////////////////////////////////////////
//! \brief The clebsch namespace is written by the authors of arXiv:1009.0437v1.
//! This code is written by the authors of that paper, but I've added a few
//! modifications to the syntax to make it consistent with my c++ style
//!
//////////////////////////////////////////////////////////////////////////////////

namespace clebsch {
    const double EPS = 1e-12;
    // binomial coefficients
    class Binomial_t {
        std::vector<int> cache;
        int N;
    public:
        int operator()(int n, int k);
    };
    extern Binomial_t binomial;
    // Eq. (19) and (25)
    class Weight {
        std::vector<int> elem;
    public:
        // the N in "SU(N)"
        const int N;
        // create a non-initialized Weight
        Weight(int N);
        // create irrep Weight of given index
        // Eq. (C2)
        Weight(int N, int index);
        // assign from another instance
        clebsch::Weight &operator=(const clebsch::Weight &w);
        // access elements of this Weight (k = 1, ..., N)
        int &operator()(int k);
        const int &operator()(int k) const;
        // compare weights
        // Eq. (C1)
        bool operator<(const Weight &w) const;
        bool operator==(const Weight &w) const;
        // element-wise sum of weights
        clebsch::Weight operator+(const Weight &w) const;
        // returns the index of this irrep Weight (index = 0, 1, ...)
        // Eq. (C2)
        int Index() const;
        // returns the dimension of this irrep Weight
        // Eq. (22)
        long long Dimension() const;
    };
    
    // Eq. (20)
    class Pattern {
        std::vector<int> elem;
    public:
        // the N in "SU(N)"
        const int N;
        // copy constructor
        Pattern(const Pattern &pat);
        // create Pattern of given index from irrep Weight
        // Eq. (C7)
        Pattern(const Weight &irrep, int index = 0);
        // access elements of this Pattern (l = 1, ..., N; k = 1, ..., l)
        int &operator()(int k, int l);
        const int &operator()(int k, int l) const;
        // find succeeding/preceding Pattern, return false if not possible
        // Eq. (C9)
        bool operator++();
        bool operator--();
        // returns the Pattern index (index = 0, ..., dimension - 1)
        // Eq. (C7)
        int Index() const;
        // returns the Pattern Weight
        // Eq. (25)
        clebsch::Weight GetWeight() const;
        // returns matrix element of lowering operator J^(l)_-
        // between this Pattern minus M^(k,l) and this Pattern
        // (l = 1, ..., N; k = 1, ..., l)
        // Eq. (28)
        double LoweringCoeff(int k, int l) const;
        // returns matrix element of raising operator J^(l)_+
        // between this Pattern plus M^(k,l) and this Pattern
        // (l = 1, ..., N; k = 1, ..., l)
        // Eq. (29)
        double RaisingCoeff(int k, int l) const;
    };
    
    class Decomposition {
        std::vector<clebsch::Weight> weights;
        std::vector<int> multiplicities;
    public:
        // the N in "SU(N)"
        const int N;
        // save given irreps for later use
        const Weight factor1, factor2;
        // construct the Decomposition of factor1 times factor2 into irreps
        // Eq. (31)
        Decomposition(const Weight &factor1, const Weight &factor2);
        // return the number of occurring irreps
        int Size() const;
        // access the occurring irreps
        // j = 0, ..., Size() - 1
        const clebsch::Weight &operator()(int j) const;
        // return the outer multiplicity of irrep in this Decomposition
        int Multiplicity(const Weight &irrep) const;
    };
    
    class IndexAdapter {
        std::vector<int> indices;
        std::vector<int> multiplicities;
    public:
        // the N in "SU(N)"
        const int N;
        // save given irreps for later use
        const int factor1, factor2;
        // construct this IndexAdapter from a given Decomposition
        IndexAdapter(const clebsch::Decomposition &decomp);
        // return the number of occurring irreps
        int Size() const;
        // access the occurring irreps
        int operator()(int j) const;
        // return the outer multiplicity of irrep in this Decomposition
        int Multiplicity(int irrep) const;
    };

    class Coefficients {
        std::map<std::vector<int>, double> clzx;
        // access Clebsch-Gordan coefficients in convenient manner
        void Set(int factor1_state,
                 int factor2_state,
                 int multiplicity_index,
                 int irrep_state,
                 double value);
        // internal functions, doing most of the work
        void HighestWeightNormalForm(); // Eq. (37)
        void ComputeHighestWeightCoeffs(); // Eq. (36)
        void ComputeLowerWeightCoeffs(int multip_index, int state, std::vector<char> &done); // Eq. (40)
    public:
        // the N in "SU(N)"
        const int N;
        // save irreps and their dimensions for later use
        const Weight factor1, factor2, irrep;
        const int factor1_dimension, factor2_dimension, irrep_dimension;
        // outer multiplicity of irrep in this Decomposition
        const int multiplicity;
        // construct all Clebsch-Gordan coefficients of this Decomposition
        Coefficients(const Weight &irrep, const Weight &factor1, const Weight &factor2);
        // access Clebsch-Gordan coefficients (read-only)
        // multiplicity_index = 0, ..., multiplicity - 1
        // factor1_state = 0, ..., factor1_dimension - 1
        // factor2_state = 0, ..., factor2_dimension - 1
        // irrep_state = 0, ..., irrep_dimension
        double operator()(int factor1_state,
                          int factor2_state,
                          int multiplicity_index,
                          int irrep_state) const;
    };
}  //  End clebsch namespace
#endif
