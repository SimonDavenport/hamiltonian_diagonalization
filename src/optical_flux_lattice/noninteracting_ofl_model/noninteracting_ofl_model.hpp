////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                 
//!                                                                             
//!	 \file
//!     This file defines a class to store the non-interacting optical
//!     flux lattice model. See e.g. PRL 109, 265301 (2012)
//!                                                                                                                                             
//!                    Copyright (C)  Simon C Davenport
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

#ifndef _NONINTERACTING_OFL_MODEL_HPP_INCLUDED_
#define _NONINTERACTING_OFL_MODEL_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../program_options/noninteracting_ofl_model_options.hpp"
#include "noninteracting_ofl_model_data.hpp"
#include "../../utilities/general/orbital_and_state_defs.hpp"
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/general/dcmplx_type_def.hpp"
#include "../../utilities/general/i_const_def.hpp"
#include "../../utilities/wrappers/lapack_wrapper.hpp"
#include "../../utilities/mathematics/lattice_vector.hpp"
#include "../../utilities/general/load_bar.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../../utilities/wrappers/io_wrapper.hpp"
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif
#include <sstream>
#include <iomanip>
#include <fstream>

namespace diagonalization
{
    namespace reciprocalLattice
    {
        //  Define the reciprocal lattice vectors used to map out a hexagonal lattice
        static const utilities::LatticeVector2D<double> G1(-0.866025403784438646763723170753,-1.5);
        static const utilities::LatticeVector2D<double> G2(0.866025403784438646763723170753,-1.5);
    };

    enum spin_t {_SPIN_UP_=0, _SPIN_DOWN_=1};            
        //!<    A data type to describe the spin index
    enum takeConjugate_t {_NO_CONJUGATE_=0, _CONJUGATE_=1};    
        //!<    A data type to flag wheter to take the complex conjugate or not

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief The NonInteractingOflModel class contains a matrix representation of
    //! the single particle optical flux lattice Hamiltonian as a 2D difference equation
    //! for Block-type eigenfunctions. The class contains funcitons to populate and
    //! diagonalize that matrix
    //!
    //! This class sets up the difference equation problem defined by taking the single
    //! particle Hamiltonian (eq. (1) in PRL 109,265301). First we perform a unitary 
    //! transformation with U = {{exp(-ik_3.r),0},{0,exp(ik_3.r)}}. This has the
    //! effect of making the Bx+iBy part (off diagonals) expressible in terms of 
    //! the reciprocal lattice vectors, at the expense of shifting the kinetic term.
    //!
    //! Then we provide the following Bloch wave function ansatz
    //!
    //! \[ |u^{n,k}_{\sigma} > = sum_{G} c^{n,k}_{G,sigma} exp(i(k+G).r) |sigma> \]
    //!
    //! Where n is the band index, k is a position in 2D k-space, {G} is the set of reciprocal
    //! lattice vectors, sigma is a spin index (up or down here) and r is the position in real
    //! space. Inserting this into Eq (1) leads to a matrix relation between the
    //! coefficients c^{n,k}_{G,sigma}. This is an infinite matrix in general,
    //! however, considered in a reduced zone scheme it becomes evident that if we only
    //! wish to determine the lowest part of the band structure then it's sufficient
    //! to solve a truncated equation with x and y cut-offs in the number of reciprocal lattice
    //! periods.
    //!
    //! The resulting coefficients are calculated for a given set of Hamiltonian
    //! parameters, and at a given position in k-space. This result can be used
    //! to build up a map of the energy bands throughout k-space, as well as the
    //! values of coefficents that can later be used to construct the interacting
    //! Hamiltonian.
    //////////////////////////////////////////////////////////////////////////////////

    class NonInteractingOflModel
    {
        private:
        friend class NonInteractingOflModelGrid;
        dcmplx* m_matrix;   //!<    Pointer to a memory address to store the matrix 
                            //!     representation of the difference equation
                            //!     relating Block wave function coefficients
        iSize_t m_dim;      //!<    Dimension of the difference equation matrix after truncation
        double* m_eigenvalues;
                            //!<    Pointer to a memory address to store a list of 
                            //!     energy eigenvalues
        dcmplx* m_blochCoefficients;
                            //!<   Pointer to a memory address to store a list of
                            //!    Bloch coefficients
        bool m_differenceMatrixPopulated;
                            //!<   Set to true once the difference equation
                            //!    matrix has been constructed
        bool m_bandsCalculated;
                            //!<   Set to true once Bloch bands have been calculated
                            //!
        NonInteractingOflModelData m_params;
                            //!<    Class internal copy of 
                            //!     the NonInteractingOflModelData struct
        utilities::LatticeVector2D<double> m_G;
                            //!<    Value in k-space where the Hamiltonian
                            //!     is to be evaluated
        int m_xBandCutOff;  //!<    Momentum space cut-off (in units of reciprocal 
                            //!     lattice vector) in the x-direction 
        int m_yBandCutOff;  //!<    Momentum space cut-off (in units of reciprocal 
                            //!     lattice vector) in the y-direction         
        iSize_t m_nbrBands;     
                            //!<    Number of bands to calculate
        //////      Private function predeclarations        ////////////////////////////
        void UpdateDifferenceDiagonal(const utilities::LatticeVector2D<double> newK);
        double EvaluateDiagonalTerm(const int l,const int m,const spin_t spinState) const;
        void MapArrayToLattice(int* x,int* y,const int i) const;
        int MapLatticeToArray(int x,int y) const;
        public: 
        NonInteractingOflModel();
        NonInteractingOflModel(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);
        NonInteractingOflModel(NonInteractingOflModelData params);
        NonInteractingOflModel(const NonInteractingOflModel &other);
        ~NonInteractingOflModel();
        void MpiSynchrnoizedCopy(const NonInteractingOflModel &other, const iSize_t nodeId, const utilities::MpiWrapper& mpi);
        void PopulateDifferenceMatrix();
        void DiagonalizeDifferenceMatrix();
	    iSize_t GetDifferenceMatrixDimension() const;
        double GetRecoilEnergy() const;
        double GetBandEnergy(const iSize_t n) const;
        void GetBlochCoefficients(dcmplx* eigenvector, const iSize_t n, const spin_t s, const takeConjugate_t conjugate) const;
        iSize_t GetCutOffX() const;
        iSize_t GetCutOffY() const;
        void PrintDifferenceMatrix();
        void MultiplyByPhase(const double phaseAngle);
        double ApplySigmaZ(const iSize_t n)  const;
        dcmplx GetSpatialWaveFunction(const spin_t spin, const double rx, const double ry, const double gridSpacing) const;
    };
}   //  End diagonalization namespace
#endif
