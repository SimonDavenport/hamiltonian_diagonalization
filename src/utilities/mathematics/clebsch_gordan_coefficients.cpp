////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 18/05/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store pre-calculated Clebsch-Gordan 
//!     coefficients. Coefficients are calculated using the algorithm described 
//!     here: arXiv:1009.0437v1. See:
//!     http://homepages.physik.uni-muenchen.de/~vondelft/Papers/ClebschGordan/
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

#include "clebsch_gordan_coefficients.hpp"

clebsch::Binomial_t clebsch::binomial;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! The following implementations are written by the authors of arXiv:1009.0437v1
//!
//////////////////////////////////////////////////////////////////////////////////

// implementation of "Binomial_t" starts here

int clebsch::Binomial_t::operator()(int n, int k) {
    if (N <= n) {
        for (cache.resize((n + 1) * (n + 2) / 2); N <= n; ++N) {
            cache[N * (N + 1) / 2] = cache[N * (N + 1) / 2 + N] = 1;
            for (int k = 1; k < N; ++k) {
                cache[N * (N + 1) / 2 + k] = cache[(N - 1) * N / 2 + k]
                                           + cache[(N - 1) * N / 2 + k - 1];
            }
        }
    }

    return cache[n * (n + 1) / 2 + k];
}

// implementation of "weight" starts here

clebsch::Weight::Weight(int N) : elem(N), N(N) {}

clebsch::Weight::Weight(int N, int index) : elem(N, 0), N(N) {
    for (int i = 0; index > 0 && i < N; ++i) {
        for (int j = 1; binomial(N - i - 1 + j, N - i - 1) <= index; j <<= 1) {
            elem[i] = j;
        }

        for (int j = elem[i] >> 1; j > 0; j >>= 1) {
            if (binomial(N - i - 1 + (elem[i] | j), N - i - 1) <= index) {
                elem[i] |= j;
            }
        }

        index -= binomial(N - i - 1 + elem[i]++, N - i - 1);
    }
}

clebsch::Weight &clebsch::Weight::operator=(const clebsch::Weight &w) {
    int &n = const_cast<int &>(N);
    elem = w.elem;
    n = w.N;
    return *this;
}

int &clebsch::Weight::operator()(int k) {
    assert(1 <= k && k <= N);
    return elem[k - 1];
}

const int &clebsch::Weight::operator()(int k) const {
    assert(1 <= k && k <= N);
    return elem[k - 1];
}

bool clebsch::Weight::operator<(const Weight &w) const {
    assert(w.N == N);
    for (int i = 0; i < N; ++i) {
        if (elem[i] - elem[N - 1] != w.elem[i] - w.elem[N - 1]) {
            return elem[i] - elem[N - 1] < w.elem[i] - w.elem[N - 1];
        }
    }
    return false;
}

bool clebsch::Weight::operator==(const Weight &w) const {
    assert(w.N == N);

    for (int i = 1; i < N; ++i) {
        if (w.elem[i] - w.elem[i - 1] != elem[i] - elem[i - 1]) {
            return false;
        }
    }

    return true;
}

clebsch::Weight clebsch::Weight::operator+(const Weight &w) const {
    Weight result(N);

    transform(elem.begin(), elem.end(), w.elem.begin(), result.elem.begin(), std::plus<int>());

    return result;
}

int clebsch::Weight::Index() const {
    int result = 0;

    for (int i = 0; elem[i] > elem[N - 1]; ++i) {
        result += binomial(N - i - 1 + elem[i] - elem[N - 1] - 1, N - i - 1);
    }

    return result;
}

long long clebsch::Weight::Dimension() const {
    long long numerator = 1, denominator = 1;

    for (int i = 1; i < N; ++i) {
        for (int j = 0; i + j < N; ++j) {
            numerator *= elem[j] - elem[i + j] + i;
            denominator *= i;
        }
    }

    return numerator / denominator;
}

// implementation of "Pattern" starts here

clebsch::Pattern::Pattern(const Pattern &p) : elem(p.elem), N(p.N) {}

clebsch::Pattern::Pattern(const Weight &irrep, int index) :
        elem((irrep.N * (irrep.N + 1)) / 2), N(irrep.N) {
    for (int i = 1; i <= N; ++i) {
        (*this)(i, N) = irrep(i);
    }

    for (int l = N - 1; l >= 1; --l) {
        for (int k = 1; k <= l; ++k) {
            (*this)(k, l) = (*this)(k + 1, l + 1);
        }
    }

    while (index-- > 0) {
        bool b = ++(*this);

        assert(b);
    }
}

int &clebsch::Pattern::operator()(int k, int l) {
    return elem[(N * (N + 1) - l * (l + 1)) / 2 + k - 1];
}

const int &clebsch::Pattern::operator()(int k, int l) const {
    return elem[(N * (N + 1) - l * (l + 1)) / 2 + k - 1];
}

bool clebsch::Pattern::operator++() {
    int k = 1, l = 1;

    while (l < N && (*this)(k, l) == (*this)(k, l + 1)) {
        if (--k == 0) {
            k = ++l;
        }
    }

    if (l == N) {
        return false;
    }

    ++(*this)(k, l);

    while (k != 1 || l != 1) {
        if (++k > l) {
            k = 1;
            --l;
        }

        (*this)(k, l) = (*this)(k + 1, l + 1);
    }

    return true;
}

bool clebsch::Pattern::operator--() {
    int k = 1, l = 1;

    while (l < N && (*this)(k, l) == (*this)(k + 1, l + 1)) {
        if (--k == 0) {
            k = ++l;
        }
    }

    if (l == N) {
        return false;
    }

    --(*this)(k, l);

    while (k != 1 || l != 1) {
        if (++k > l) {
            k = 1;
            --l;
        }

        (*this)(k, l) = (*this)(k, l + 1);
    }

    return true;
}

int clebsch::Pattern::Index() const {
    int result = 0;

    for (Pattern p(*this); --p; ++result) {}

    return result;
}

clebsch::Weight clebsch::Pattern::GetWeight() const {
    clebsch::Weight result(N);

    for (int prev = 0, l = 1; l <= N; ++l) {
        int now = 0;

        for (int k = 1; k <= l; ++k) {
            now += (*this)(k, l);
        }

        result(l) = now - prev;
        prev = now;
    }

    return result;
}

double clebsch::Pattern::LoweringCoeff(int k, int l) const {
    double result = 1.0;

    for (int i = 1; i <= l + 1; ++i) {
        result *= (*this)(i, l + 1) - (*this)(k, l) + k - i + 1;
    }
    
    for (int i = 1; i <= l - 1; ++i) {
        result *= (*this)(i, l - 1) - (*this)(k, l) + k - i;
    }

    for (int i = 1; i <= l; ++i) {
        if (i == k) continue;
        result /= (*this)(i, l) - (*this)(k, l) + k - i + 1;
        result /= (*this)(i, l) - (*this)(k, l) + k - i;
    }

    return std::sqrt(-result);
}

double clebsch::Pattern::RaisingCoeff(int k, int l) const {
    double result = 1.0;

    for (int i = 1; i <= l + 1; ++i) {
        result *= (*this)(i, l + 1) - (*this)(k, l) + k - i;
    }

    for (int i = 1; i <= l - 1; ++i) {
        result *= (*this)(i, l - 1) - (*this)(k, l) + k - i - 1;
    }

    for (int i = 1; i <= l; ++i) {
        if (i == k) continue;
        result /= (*this)(i, l) - (*this)(k, l) + k - i;
        result /= (*this)(i, l) - (*this)(k, l) + k - i  - 1;
    }

    return std::sqrt(-result);
}

// implementation of "Decomposition" starts here

clebsch::Decomposition::Decomposition(const Weight &factor1, const Weight &factor2) :
        N(factor1.N), factor1(factor1), factor2(factor2) {
    assert(factor1.N == factor2.N);
    std::vector<clebsch::Weight> result;
    Pattern low(factor1), high(factor1);
    Weight trial(factor2);
    int k = 1, l = N;

    do {
        while (k <= N) {
            --l;
            if (k <= l) {
                low(k, l) = std::max(high(k + N - l, N), high(k, l + 1) + trial(l + 1) - trial(l));
                high(k, l) = high(k, l + 1);
                if (k > 1 && high(k, l) > high(k - 1, l - 1)) {
                    high(k, l) = high(k - 1, l - 1);
                }
                if (l > 1 && k == l && high(k, l) > trial(l - 1) - trial(l)) {
                    high(k, l) = trial(l - 1) - trial(l);
                }
                if (low(k, l) > high(k, l)) {
                    break;
                }
                trial(l + 1) += high(k, l + 1) - high(k, l);
            } else {
                trial(l + 1) += high(k, l + 1);
                ++k;
                l = N;
            }
        }

        if (k > N) {
            result.push_back(trial);
            for (int i = 1; i <= N; ++i) {
                result.back()(i) -= result.back()(N);
            }
        } else {
            ++l;
        }

        while (k != 1 || l != N) {
            if (l == N) {
                l = --k - 1;
                trial(l + 1) -= high(k, l + 1);
            } else if (low(k, l) < high(k, l)) {
                --high(k, l);
                ++trial(l + 1);
                break;
            } else {
                trial(l + 1) -= high(k, l + 1) - high(k, l);
            }
            ++l;
        }
    } while (k != 1 || l != N);

    sort(result.begin(), result.end());
    for (std::vector<clebsch::Weight>::iterator it = result.begin(); it != result.end(); ++it) {
        if (it != result.begin() && *it == weights.back()) {
            ++multiplicities.back();
        } else {
            weights.push_back(*it);
            multiplicities.push_back(1);
        }
    }
}

int clebsch::Decomposition::Size() const {
    return weights.size();
}

const clebsch::Weight &clebsch::Decomposition::operator()(int j) const {
    return weights[j];
}

int clebsch::Decomposition::Multiplicity(const Weight &irrep) const {
    assert(irrep.N == N);
    std::vector<clebsch::Weight>::const_iterator it
        = std::lower_bound(weights.begin(), weights.end(), irrep);

    return it != weights.end() && *it == irrep ? multiplicities[it - weights.begin()] : 0;
}

// implementation of "IndexAdapter" starts here

clebsch::IndexAdapter::IndexAdapter(const clebsch::Decomposition &decomp) :
        N(decomp.N),
        factor1(decomp.factor1.Index()),
        factor2(decomp.factor2.Index()) {
    for (int i = 0, s = decomp.Size(); i < s; ++i) {
        indices.push_back(decomp(i).Index());
        multiplicities.push_back(decomp.Multiplicity(decomp(i)));
    }
}

int clebsch::IndexAdapter::Size() const {
    return indices.size();
}

int clebsch::IndexAdapter::operator()(int j) const {
    return indices[j];
}

int clebsch::IndexAdapter::Multiplicity(int irrep) const {
    std::vector<int>::const_iterator it = std::lower_bound(indices.begin(), indices.end(), irrep);

    return it != indices.end() && *it == irrep ? multiplicities[it - indices.begin()] : 0;
}

// implementation of "Coefficients" starts here

void clebsch::Coefficients::Set(int factor1_state,
                                int factor2_state,
                                int multiplicity_index,
                                int irrep_state,
                                double value) {
    assert(0 <= factor1_state && factor1_state < factor1_dimension);
    assert(0 <= factor2_state && factor2_state < factor2_dimension);
    assert(0 <= multiplicity_index && multiplicity_index < multiplicity);
    assert(0 <= irrep_state && irrep_state < irrep_dimension);

    int coefficient_label[] = { factor1_state,
                                factor2_state,
                                multiplicity_index,
                                irrep_state };
    clzx[std::vector<int>(coefficient_label, coefficient_label
            + sizeof coefficient_label / sizeof coefficient_label[0])] = value;
}

void clebsch::Coefficients::HighestWeightNormalForm() {
    int hws = irrep_dimension - 1;

    // bring CGCs into reduced row echelon form
    for (int h = 0, i = 0; h < multiplicity - 1 && i < factor1_dimension; ++i) {
        for (int j = 0; h < multiplicity - 1 && j < factor2_dimension; ++j) {
            int k0 = h;

            for (int k = h + 1; k < multiplicity; ++k) {
                if (fabs((*this)(i, j, k, hws)) > fabs((*this)(i, j, k0, hws))) {
                    k0 = k;
                }
            }

            if ((*this)(i, j, k0, hws) < -EPS) {
                for (int i2 = i; i2 < factor1_dimension; ++i2) {
                    for (int j2 = i2 == i ? j : 0; j2 < factor2_dimension; ++j2) {
                        Set(i2, j2, k0, hws, -(*this)(i2, j2, k0, hws));
                    }
                }
            } else if ((*this)(i, j, k0, hws) < EPS) {
                continue;
            }

            if (k0 != h) {
                for (int i2 = i; i2 < factor1_dimension; ++i2) {
                    for (int j2 = i2 == i ? j : 0; j2 < factor2_dimension; ++j2) {
                        double x = (*this)(i2, j2, k0, hws);
                        Set(i2, j2, k0, hws, (*this)(i2, j2, h, hws));
                        Set(i2, j2, h, hws, x);
                    }
                }
            }

            for (int k = h + 1; k < multiplicity; ++k) {
                for (int i2 = i; i2 < factor1_dimension; ++i2) {
                    for (int j2 = i2 == i ? j : 0; j2 < factor2_dimension; ++j2) {
                        Set(i2, j2, k, hws, (*this)(i2, j2, k, hws) - (*this)(i2, j2, h, hws)
                                * (*this)(i, j, k, hws) / (*this)(i, j, h, hws));
                    }
                }
            }

            // next 3 lines not strictly necessary, might improve numerical stability
            for (int k = h + 1; k < multiplicity; ++k) {
                Set(i, j, k, hws, 0.0);
            }

            ++h;
        }
    }

    // Gram-Schmidt orthonormalization
    for (int h = 0; h < multiplicity; ++h) {
        for (int k = 0; k < h; ++k) {
            double overlap = 0.0;
            for (int i = 0; i < factor1_dimension; ++i) {
                for (int j = 0; j < factor2_dimension; ++j) {
                    overlap += (*this)(i, j, h, hws) * (*this)(i, j, k, hws);
                }
            }

            for (int i = 0; i < factor1_dimension; ++i) {
                for (int j = 0; j < factor2_dimension; ++j) {
                    Set(i, j, h, hws, (*this)(i, j, h, hws) - overlap * (*this)(i, j, k, hws));
                }
            }
        }

        double norm = 0.0;
        for (int i = 0; i < factor1_dimension; ++i) {
            for (int j = 0; j < factor2_dimension; ++j) {
                norm += (*this)(i, j, h, hws) * (*this)(i, j, h, hws);
            }
        }
        norm = std::sqrt(norm);

        for (int i = 0; i < factor1_dimension; ++i) {
            for (int j = 0; j < factor2_dimension; ++j) {
                Set(i, j, h, hws, (*this)(i, j, h, hws) / norm);
            }
        }
    }
}

void clebsch::Coefficients::ComputeHighestWeightCoeffs() {
    if (multiplicity == 0) {
        return;
    }

    std::vector<std::vector<int> > map_coeff(factor1_dimension,
                                             std::vector<int>(factor2_dimension, -1));
    std::vector<std::vector<int> > map_states(factor1_dimension,
                                              std::vector<int>(factor2_dimension, -1));
    int n_coeff = 0, n_states = 0;
    Pattern p(factor1, 0);

    for (int i = 0; i < factor1_dimension; ++i, ++p) {
        Weight pw(p.GetWeight());
        Pattern q(factor2, 0);
        for (int j = 0; j < factor2_dimension; ++j, ++q) {
            if (pw + q.GetWeight() == irrep) {
                map_coeff[i][j] = n_coeff++;
            }
        }
    }

    if (n_coeff == 1) {
        for (int i = 0; i < factor1_dimension; ++i) {
            for (int j = 0; j < factor2_dimension; ++j) {
                if (map_coeff[i][j] >= 0) {
                    Set(i, j, 0, irrep_dimension - 1, 1.0);
                    return;
                }
            }
        }
    }

    double *hw_system = new double[n_coeff * (factor1_dimension * factor2_dimension)];
    assert(hw_system != NULL);
    memset(hw_system, 0, n_coeff * (factor1_dimension * factor2_dimension) * sizeof (double));

    Pattern r(factor1, 0);
    for (int i = 0; i < factor1_dimension; ++i, ++r) {
        Pattern q(factor2, 0);

        for (int j = 0; j < factor2_dimension; ++j, ++q) {
            if (map_coeff[i][j] >= 0) {
                for (int l = 1; l <= N - 1; ++l) {
                    for (int k = 1; k <= l; ++k) {
                        if ((k == 1 || r(k, l) + 1 <= r(k - 1, l - 1)) && r(k, l) + 1 <= r(k, l + 1)) {
                            ++r(k, l);
                            int h = r.Index();
                            --r(k, l);

                            if (map_states[h][j] < 0) {
                                map_states[h][j] = n_states++;
                            }

                            hw_system[n_coeff * map_states[h][j] + map_coeff[i][j]]
                                += r.RaisingCoeff(k, l);
                        }

                        if ((k == 1 || q(k, l) + 1 <= q(k - 1, l - 1)) && q(k, l) + 1 <= q(k, l + 1)) {
                            ++q(k, l);
                            int h = q.Index();
                            --q(k, l);

                            if (map_states[i][h] < 0) {
                                map_states[i][h] = n_states++;
                            }


                            hw_system[n_coeff * map_states[i][h] + map_coeff[i][j]]
                                += q.RaisingCoeff(k, l);
                        }
                    }
                }
            }
        }
    }

    int lwork = -1, info;
    double worksize;

    double *singval = new double[std::min(n_coeff, n_states)];
    assert(singval != NULL);
    double *singvec = new double[n_coeff * n_coeff];
    assert(singvec != NULL);

    dgesvd_("A",
            "N",
            &n_coeff,
            &n_states,
            hw_system,
            &n_coeff,
            singval,
            singvec,
            &n_coeff,
            NULL,
            &n_states,
            &worksize,
            &lwork,
            &info);
    assert(info == 0);

    lwork = worksize;
    double *work = new double[lwork];
    assert(work != NULL);

    dgesvd_("A",
            "N",
            &n_coeff,
            &n_states,
            hw_system,
            &n_coeff,
            singval,
            singvec,
            &n_coeff,
            NULL,
            &n_states,
            work,
            &lwork,
            &info);
    assert(info == 0);

    for (int i = 0; i < multiplicity; ++i) {
        for (int j = 0; j < factor1_dimension; ++j) {
            for (int k = 0; k < factor2_dimension; ++k) {
                if (map_coeff[j][k] >= 0) {
                    double x = singvec[n_coeff * (n_coeff - 1 - i) + map_coeff[j][k]];

                    if (fabs(x) > EPS) {
                        Set(j, k, i, irrep_dimension - 1, x);
                    }
                }
            }
        }
    }

    // uncomment next line to bring highest-weight coefficients into "normal form"
    // HighestWeightNormalForm();

    delete[] work;
    delete[] singvec;
    delete[] singval;
    delete[] hw_system;
}

void clebsch::Coefficients::ComputeLowerWeightCoeffs(int multip_index,
                                                        int state,
                                                        std::vector<char> &done) {
    Weight statew(Pattern(irrep, state).GetWeight());
    Pattern p(irrep, 0);
    std::vector<int> map_parent(irrep_dimension, -1),
                     map_multi(irrep_dimension, -1),
                     which_l(irrep_dimension, -1);
    int n_parent = 0, n_multi = 0;

    for (int i = 0; i < irrep_dimension; ++i, ++p) {
        Weight v(p.GetWeight());

        if (v == statew) {
            map_multi[i] = n_multi++;
        } else for (int l = 1; l < N; ++l) {
            --v(l);
            ++v(l + 1);
            if (v == statew) {
                map_parent[i] = n_parent++;
                which_l[i] = l;
                if (!done[i]) {
                    ComputeLowerWeightCoeffs(multip_index, i, done);
                }
                break;
            }
            --v(l + 1);
            ++v(l);
        }
    }

    double *irrep_coeffs = new double[n_parent * n_multi];
    assert(irrep_coeffs != NULL);
    memset(irrep_coeffs, 0, n_parent * n_multi * sizeof (double));

    double *prod_coeffs = new double[n_parent * factor1_dimension * factor2_dimension];
    assert(prod_coeffs != NULL);
    memset(prod_coeffs, 0, n_parent * factor1_dimension * factor2_dimension * sizeof (double));

    std::vector<std::vector<int> > map_prodstat(factor1_dimension,
                                                std::vector<int>(factor2_dimension, -1));
    int n_prodstat = 0;

    Pattern r(irrep, 0);
    for (int i = 0; i < irrep_dimension; ++i, ++r) {
        if (map_parent[i] >= 0) {
            for (int k = 1, l = which_l[i]; k <= l; ++k) {
                if (r(k, l) > r(k + 1, l + 1) && (k == l || r(k, l) > r(k, l - 1))) {
                    --r(k, l);
                    int h = r.Index();
                    ++r(k, l);

                    irrep_coeffs[n_parent * map_multi[h] + map_parent[i]] += r.LoweringCoeff(k, l);
                }
            }

            Pattern q1(factor1, 0);
            for (int j1 = 0; j1 < factor1_dimension; ++j1, ++q1) {
                Pattern q2(factor2, 0);

                for (int j2 = 0; j2 < factor2_dimension; ++j2, ++q2) {
                    if (std::fabs((*this)(j1, j2, multip_index, i)) > EPS) {
                        for (int k = 1, l = which_l[i]; k <= l; ++k) {
                            if (q1(k, l) > q1(k + 1, l + 1) && (k == l || q1(k, l) > q1(k, l - 1))) {
                                --q1(k, l);
                                int h = q1.Index();
                                ++q1(k, l);

                                if (map_prodstat[h][j2] < 0) {
                                    map_prodstat[h][j2] = n_prodstat++;
                                }

                                prod_coeffs[n_parent * map_prodstat[h][j2] + map_parent[i]] +=
                                        (*this)(j1, j2, multip_index, i) * q1.LoweringCoeff(k, l);
                            }

                            if (q2(k, l) > q2(k + 1, l + 1) && (k == l || q2(k, l) > q2(k, l - 1))) {
                                --q2(k, l);
                                int h = q2.Index();
                                ++q2(k, l);

                                if (map_prodstat[j1][h] < 0) {
                                    map_prodstat[j1][h] = n_prodstat++;
                                }

                                prod_coeffs[n_parent * map_prodstat[j1][h] + map_parent[i]] +=
                                        (*this)(j1, j2, multip_index, i) * q2.LoweringCoeff(k, l);
                            }
                        }
                    }
                }
            }
        }
    }

    double worksize;
    int lwork = -1, info;

    dgels_("N",
           &n_parent,
           &n_multi,
           &n_prodstat,
           irrep_coeffs,
           &n_parent,
           prod_coeffs,
           &n_parent,
           &worksize,
           &lwork,
           &info);
    assert(info == 0);

    lwork = worksize;
    double *work = new double[lwork];
    assert(work != NULL);

    dgels_("N",
           &n_parent,
           &n_multi,
           &n_prodstat,
           irrep_coeffs,
           &n_parent,
           prod_coeffs,
           &n_parent,
           work,
           &lwork,
           &info);
    assert(info == 0);

    for (int i = 0; i < irrep_dimension; ++i) {
        if (map_multi[i] >= 0) {
            for (int j = 0; j < factor1_dimension; ++j) {
                for (int k = 0; k < factor2_dimension; ++k) {
                    if (map_prodstat[j][k] >= 0) {
                        double x = prod_coeffs[n_parent * map_prodstat[j][k] + map_multi[i]];

                        if (fabs(x) > EPS) {
                            Set(j, k, multip_index, i, x);
                        }
                    }
                }
            }

            done[i] = true;
        }
    }

    delete[] work;
    delete[] prod_coeffs;
    delete[] irrep_coeffs;
}

clebsch::Coefficients::Coefficients(const Weight &irrep, const Weight &factor1, const Weight &factor2) :
        N(irrep.N),
        factor1(factor1),
        factor2(factor2),
        irrep(irrep),
        factor1_dimension(factor1.Dimension()),
        factor2_dimension(factor2.Dimension()),
        irrep_dimension(irrep.Dimension()),
        multiplicity(Decomposition(factor1, factor2).Multiplicity(irrep)) {
    assert(factor1.N == irrep.N);
    assert(factor2.N == irrep.N);

    ComputeHighestWeightCoeffs();

    for (int i = 0; i < multiplicity; ++i) {
        std::vector<char> done(irrep_dimension, 0);
        done[irrep_dimension - 1] = true;
        for (int j = irrep_dimension - 1; j >= 0; --j) {
            if (!done[j]) {
                ComputeLowerWeightCoeffs(i, j, done);
            }
        }
    }
}

double clebsch::Coefficients::operator()(int factor1_state,
                                         int factor2_state,
                                         int multiplicity_index,
                                         int irrep_state) const {
    assert(0 <= factor1_state && factor1_state < factor1_dimension);
    assert(0 <= factor2_state && factor2_state < factor2_dimension);
    assert(0 <= multiplicity_index && multiplicity_index < multiplicity);
    assert(0 <= irrep_state && irrep_state < irrep_dimension);

    int coefficient_label[] = { factor1_state,
                                factor2_state,
                                multiplicity_index,
                                irrep_state };
    std::map<std::vector<int>, double>::const_iterator it(
            clzx.find(std::vector<int>(coefficient_label, coefficient_label
                    + sizeof coefficient_label / sizeof coefficient_label[0])));

    return it != clzx.end() ? it->second : 0.0;
}
