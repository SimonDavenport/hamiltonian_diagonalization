////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                
//!                                                                             
//!	 \file
//!     This file defines a class to store the non-interacting optical
//!     flux lattice model. See e.g. PRL 109, 265301 (2012)	 
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "noninteracting_ofl_model.hpp"

namespace diagonalization
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Default constructor for the NonInteractingOflModel class
    //! 
    //! Initializes the class internal pointers to NULL and initializes
    //! the class bools to false
    ////////////////////////////////////////////////////////////////////////////////
    NonInteractingOflModel::NonInteractingOflModel()
        :   m_matrix(0),
            m_eigenvalues(0),
            m_blochCoefficients(0),
            m_differenceMatrixPopulated(false),
            m_bandsCalculated(false),
            m_xBandCutOff(6),
            m_yBandCutOff(6),
            m_nbrBands(2)
    {
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor from command line arguments
    //! 
    //! Also initializes the class internal pointers to NULL and initializes
    //! the class bools to false
    ////////////////////////////////////////////////////////////////////////////////
    NonInteractingOflModel::NonInteractingOflModel(	
        boost::program_options::variables_map* optionList,  
                                            //!<    Parsed command line argument list
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
        :   m_matrix(0),
            m_eigenvalues(0),
            m_blochCoefficients(0),
            m_differenceMatrixPopulated(false),
            m_bandsCalculated(false)
    {
        m_params = NonInteractingOflModelData(optionList,mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            m_xBandCutOff = (*optionList)["x-cut"].as<iSize_t>();
            m_yBandCutOff = (*optionList)["y-cut"].as<iSize_t>(); 
            m_nbrBands    = (*optionList)["nbr-bands"].as<iSize_t>();
        }
        mpi.Sync(&m_xBandCutOff, 1, 0);
        mpi.Sync(&m_yBandCutOff, 1, 0);
        mpi.Sync(&m_nbrBands, 1, 0);
	    return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor for the NonInteractingOflModel class
    //! 
    //! Uses a pre-declared SingleParticleParameters struct to initialize the class.
    ////////////////////////////////////////////////////////////////////////////////
    NonInteractingOflModel::NonInteractingOflModel(
        NonInteractingOflModelData params)   //!<    An instance of the SingleParticleParameters 
                                             //!     struct to initialize the class
    :   m_matrix(0),
        m_eigenvalues(0),
        m_blochCoefficients(0),
        m_differenceMatrixPopulated(false),
        m_bandsCalculated(false),
        m_params(params),
        m_xBandCutOff(6),
        m_yBandCutOff(6),
        m_nbrBands(2)
    {}
    
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Copy constructor for the NonInteractingOflModel class
    //! 
    //! Copies all the arrays that the class contains, if they have been populated
    ////////////////////////////////////////////////////////////////////////////////
    NonInteractingOflModel::NonInteractingOflModel(
        const NonInteractingOflModel &other) //!<    Address of the class being copied
        :   m_matrix(0),
            m_dim(other.m_dim),
            m_eigenvalues(0),
            m_blochCoefficients(0),
            m_differenceMatrixPopulated(other.m_differenceMatrixPopulated),
            m_bandsCalculated(other.m_bandsCalculated),
            m_params(other.m_params),
            m_G(other.m_G),
            m_xBandCutOff(other.m_xBandCutOff),
            m_yBandCutOff(other.m_yBandCutOff),
            m_nbrBands(other.m_nbrBands)
    {
        if(m_differenceMatrixPopulated)
        {
            m_matrix = new (std::nothrow) dcmplx[m_dim*m_dim];
            memcpy (m_matrix, other.m_matrix, sizeof(dcmplx)*m_dim*m_dim);
        }
        if(m_bandsCalculated)
        {
            //  Allocate memory to store eigenvectors for lowest the specified number of bands
            m_blochCoefficients = new (std::nothrow) dcmplx[m_dim*m_nbrBands];
            //  Allocate memory to store eigenvalues for lowest bands
            m_eigenvalues = new (std::nothrow) double[m_nbrBands];
            memcpy(m_blochCoefficients, other.m_blochCoefficients, sizeof(dcmplx)*m_dim*m_nbrBands);
            memcpy(m_eigenvalues, other.m_eigenvalues, sizeof(double)*m_nbrBands);
        }
        return;
    }

    //!
    //! \brief Copy given class and synchronize over all nodes  
    //!
    void NonInteractingOflModel::MpiSynchrnoizedCopy(
        const NonInteractingOflModel &other, //!<    Object address to be copied
        const iSize_t nodeId,   //!<    Node id to sync with
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        m_params = other.m_params;
        m_dim = other.m_dim;
        m_xBandCutOff = other.m_xBandCutOff;
        m_yBandCutOff = other.m_yBandCutOff;
        m_nbrBands    = other.m_nbrBands;
        m_differenceMatrixPopulated = other.m_differenceMatrixPopulated;
        m_bandsCalculated = other.m_bandsCalculated;
        m_G = other.m_G;
        m_params.MpiSynchronize(nodeId, mpi);
        MPI_Bcast(&m_dim, 1, mpi.GetType<iSize_t>(), nodeId, mpi.m_comm);
        MPI_Bcast(&m_xBandCutOff, 1, mpi.GetType<int>(), nodeId, mpi.m_comm);
        MPI_Bcast(&m_yBandCutOff, 1, mpi.GetType<int>(), nodeId, mpi.m_comm);
        MPI_Bcast(&m_nbrBands, 1, mpi.GetType<iSize_t>(), nodeId, mpi.m_comm);
        MPI_Bcast(&m_differenceMatrixPopulated, 1, mpi.GetType<bool>(), nodeId, mpi.m_comm);
        MPI_Bcast(&m_bandsCalculated, 1, mpi.GetType<bool>(), nodeId, mpi.m_comm);
        m_G.MpiSync(nodeId, mpi);
        if(m_differenceMatrixPopulated)
        {
            if(0 == nodeId)
            {
                std::cerr<<"ERROR: NonInteractingOflModel::MpiSynchrnoizedCopy:"
                <<" MATRIX COPY NOT IMPLEMENTED"<<std::endl;
            }
        }
        if(m_bandsCalculated)
        {
            //  Allocate memory to store eigenvectors for lowest the specified number of bands
            m_blochCoefficients = new (std::nothrow) dcmplx[m_dim*m_nbrBands];
            //  Allocate memory to store eigenvalues for lowest bands
            m_eigenvalues = new (std::nothrow) double[m_nbrBands];
            //  Direct copy on the current node
	        if(mpi.m_id == (int)nodeId)
	        {
                memcpy(m_blochCoefficients, other.m_blochCoefficients, sizeof(dcmplx)*m_dim*m_nbrBands);
                memcpy(m_eigenvalues, other.m_eigenvalues, sizeof(double)*m_nbrBands);
            }
	        MPI_Barrier(mpi.m_comm);
            MPI_Bcast(m_blochCoefficients, m_dim*m_nbrBands, mpi.GetType<dcmplx>(), nodeId, mpi.m_comm);
            MPI_Bcast(m_eigenvalues, m_nbrBands, MPI_DOUBLE, nodeId, mpi.m_comm); 
        }
        return;
    }

    //!
    //! Destructor for the NonInteractingOflModel class
    //! 
    NonInteractingOflModel::~NonInteractingOflModel()
    {
        if(m_matrix!=0)              delete[] m_matrix;
        if(m_eigenvalues!=0)         delete[] m_eigenvalues;
        if(m_blochCoefficients!=0)   delete[] m_blochCoefficients;
        return;
    }

    //!
    //! Display the value of the difference equation matrix on screen
    //! 
    void NonInteractingOflModel::PrintDifferenceMatrix()
    {
        if(m_differenceMatrixPopulated)
        {    
            utilities::cout.SecondaryOutput().precision(6);
            const int lineWidth=18;    
            utilities::cout.SecondaryOutput()<<std::endl;
            if(m_dim<10)    //  for small matrices
            {
                for(iSize_t i=0; i<m_dim; ++i)
                {
	                utilities::cout.SecondaryOutput()<<"\t";
	                for(iSize_t j=0; j<m_dim; ++j)
	                {
		                utilities::cout.SecondaryOutput()<<std::setw(lineWidth)<<m_matrix[i*m_dim+j]<<" ";
	                    if((j+1)%m_dim==0)
	                    {
	                        utilities::cout.SecondaryOutput()<<"|";
	                    }
	                }
	                utilities::cout.SecondaryOutput()<<std::endl;
	            }
	        }
	        else
	        {
	            //  only print non-zeros if the matrix is large
	            for(iSize_t i=0; i<m_dim; ++i)
                {
	                for(iSize_t j=0; j<m_dim; ++j)
	                {
	                    if(abs(m_matrix[i*m_dim+j])>0.00000001)
	                    {
	                        utilities::cout.SecondaryOutput()<<"\t m["<<i<<","<<j<<"] = "<<m_matrix[i*m_dim+j]<<std::endl;
	                    }   
	                }
	            }
	        }
	    }
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Populate the difference equation relating Bloch coefficients
    //! with different reciprocal lattice indices
    //!
    //!	The matrix is stored with spin up /spin down components alternating
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModel::PopulateDifferenceMatrix()
    {
        if(m_differenceMatrixPopulated)
        {
		    delete[] m_matrix;
        }
        //  The difference equation includes contributions form both positive and 
        //  negative reciprocal lattice points.
        //  The dimension of the difference equation matrix is given by the number of
	    //	points in 2D k-space multiplied by the number of spin species
        m_dim = this->GetDifferenceMatrixDimension();
	    m_matrix = new (std::nothrow) dcmplx[m_dim*m_dim];
        const dcmplx alpha = m_params.m_V0*(3*cos(m_params.m_theta)*cos(m_params.m_theta)-1.0)/2.0;
        const dcmplx beta  = dcmplx(0,-sqrt(3.0)*m_params.m_V0*sin(m_params.m_theta)*sin(m_params.m_theta)/2.0);
        const dcmplx gamma = m_params.m_epsilon*m_params.m_V0*cos(m_params.m_theta);
        const dcmplx alphaPlusBeta = alpha + beta;
        const dcmplx alphaMinusBeta = alpha - beta;
        //  Populate the matrix (upper triangular elements only)
        for(iSize_t i=0; i<m_dim; ++i)
        {
            int lRow,mRow;   
            this->MapArrayToLattice(&lRow, &mRow, i/2);    
            for(iSize_t j=i;j<m_dim;++j)
            {
                //  Keep track of the book keeping in assigning zero /non-zero 
                //  matrix elements
                int lCol,mCol;
                this->MapArrayToLattice(&lCol, &mCol, j/2);
                const int deltaL = lCol - lRow;
                const int deltaM = mCol - mRow; 
                const spin_t spinStateI = i&1 ? _SPIN_DOWN_ : _SPIN_UP_;
                const spin_t spinStateJ = j&1 ? _SPIN_DOWN_ : _SPIN_UP_;
                dcmplx element = 0.0;
                if(i==j)
                {            
                    element = this->EvaluateDiagonalTerm(lCol,mCol,spinStateI);  
                }
                else if(std::abs(deltaL) > 1 || std::abs(deltaM) > 1)
                {
                    element = 0.0;
                }
                else
                {        
			        if(spinStateI==_SPIN_UP_ && spinStateJ==_SPIN_UP_)
			        {
				        //	spin up - spin up terms
				        if(deltaL == 0 && deltaM == 1)
                        {
                            element = alphaPlusBeta;
                        }
                        else if(deltaL == 0 && deltaM == -1)
                        {
                            element = alphaMinusBeta;
                        }
                        else if(deltaL == 1 && deltaM == 0)
                        {
                            element = alphaMinusBeta;
                        }
                        else if(deltaL == -1 && deltaM == 0)
                        {
                            element = alphaPlusBeta;
                        }
                        else if(deltaL == 1 && deltaM == -1)
                        {
                            element = alphaPlusBeta;
                        }
                        else if(deltaL == -1 && deltaM == 1)
                        {
                            element = alphaMinusBeta;
                        }
			        }
			        else if(spinStateI==_SPIN_UP_ && spinStateJ==_SPIN_DOWN_)
			        {
				        //	spin up - spin down terms
				        if(deltaL == 0 && deltaM == -1)
                        {
                            element = gamma;
                        }
                        else if(deltaL == -1 && deltaM == 0)
                        {
                            element = gamma;
                        }
                        else if(deltaL == 0 && deltaM == 0)
                        {
                            element = gamma;
                        }
			        }
			        else if(spinStateI==_SPIN_DOWN_ && spinStateJ==_SPIN_UP_)
			        {
				        //	spin down - spin up terms
				        if(deltaL == 0 && deltaM == 1)
                        {
                            element = gamma;
                        }
                        else if(deltaL == 1 && deltaM == 0)
                        {
                            element = gamma;
                        }
                        else if(deltaL == 0 && deltaM == 0)
                        {
                            element = gamma;
                        }
			        }
			        else
			        {
				        //	spin down spin down terms
				        if(deltaL == 0 && deltaM == 1)
                        {
                            element = alphaMinusBeta;
                        }
                        else if(deltaL == 0 && deltaM == -1)
                        {
                            element = alphaPlusBeta;
                        }
                        else if(deltaL == 1 && deltaM == 0)
                        {
                            element = alphaPlusBeta;
                        }
                        else if(deltaL == -1 && deltaM == 0)
                        {
                            element = alphaMinusBeta;
                        }
                        else if(deltaL == 1 && deltaM == -1)
                        {
                            element = alphaMinusBeta;
                        }
                        else if(deltaL == -1 && deltaM == 1)
                        {
                            element = alphaPlusBeta;
                        }
			        } 
                }
	            m_matrix[i*m_dim+j] = element;
            }
        }
        m_differenceMatrixPopulated = true;
    }

    //!
    //! If we change only the point in k-space then we don't need to recalculate 
    //! the off-diagonal elements of the difference matrix. 
    //!
    void NonInteractingOflModel::UpdateDifferenceDiagonal(
        const utilities::LatticeVector2D<double> newK)       //!<    Updated lattice vector
    {
        if(m_differenceMatrixPopulated)
        {
            //  Update the class internal G
            m_G = newK;
            //  Recalculate the diagonal element of the difference matrix. 
            for(iSize_t i=0; i<m_dim; ++i)
            {
                int lCol,mCol;
                const spin_t spinStateI = i&1 ? _SPIN_DOWN_ : _SPIN_UP_;
                this->MapArrayToLattice(&lCol, &mCol, i/2);
                const double element = this->EvaluateDiagonalTerm(lCol, mCol, spinStateI);
	            m_matrix[i*m_dim+i] = element;
            }
        }
        else
        {}
        return;
    }

    //!
    //! Evaluate diagonal term in difference matrix
    //!
    double NonInteractingOflModel::EvaluateDiagonalTerm(
	    const int l,			//!<	number of G1 lattice vectors
	    const int m,			//!<	number of G2 lattice vectors
	    const spin_t spinState)	//!<	Spin state of the matrix element
	    const 
    {
	    //  Diagonal terms only contain kinetic energy 
	    //  (shifted by - 0.5 Sigma_z due to the gauge transformation)
	    double temp1 = m_G.GetKx() + l*reciprocalLattice::G1.GetKx() + m*reciprocalLattice::G2.GetKx();
	    double temp2 = m_G.GetKy() + l*reciprocalLattice::G1.GetKy() + m*reciprocalLattice::G2.GetKy();
	    if(spinState==_SPIN_UP_)
	    {
		     temp2 += 1.0/2.0;
	    }
	    else
	    {
		     temp2 -= 1.0/2.0;
	    }
	    return this->GetRecoilEnergy()*(temp1*temp1 + temp2*temp2 + 6.0*m_params.m_V0);
    }

    //!
    //! Diagonalize the difference equation to determine the energy bands
    //! and the values of the Bloch coefficients
    //!
    void NonInteractingOflModel::DiagonalizeDifferenceMatrix()
    {
        if(m_differenceMatrixPopulated)
        {
            //  Remove any previous memory allocation
            if(m_blochCoefficients!=0)  delete[] m_blochCoefficients;
            if(m_eigenvalues!=0)        delete[] m_eigenvalues;
            //  Allocate memory to store eigenvectors for lowest the specified number of bands
            m_blochCoefficients = new (std::nothrow) dcmplx[m_dim*m_nbrBands];
            //  Allocate memory to store eigenvalues for lowest bands           
            m_eigenvalues = new (std::nothrow) double[m_nbrBands];
            if(m_params.m_V0==0.0)
            {}
            else
            {   
	            //  Generate a copy of the matrix to be diagonalized
	            //  (since the diagonalization routine will destroy the 
	            //  original matrix, which we want to keep)
	            dcmplx *matrix = new dcmplx[m_dim*m_dim];
	            memcpy(matrix, m_matrix, m_dim*m_dim*sizeof(dcmplx));
	            utilities::linearAlgebra::DiagonalizeSymmetricMatrix<dcmplx>(matrix, m_blochCoefficients, m_eigenvalues, 
	                                                                         m_dim, m_nbrBands, 'U');
	            delete[] matrix;
            }
        }
        else
        {}
        m_bandsCalculated = true;
        return;
    }

    //!
    //! Return the nth eigenvalue obtained from the difference equation
    //!
    double NonInteractingOflModel::GetBandEnergy(
        const iSize_t n)        //!<    The energy band index requested
        const
    {
        if(n>=0 && n<m_nbrBands)
        {
            if(m_bandsCalculated)
            {
                return m_eigenvalues[n];
            }
        }    
        return 0;
    }
                     
    //!
    //!  Map an array index i to a pair of reciprocal lattice indices x,y
    //!  such that x,y span a reciprocal lattice in both +ve and -ve directions
    //!
    void NonInteractingOflModel::MapArrayToLattice(
        int* x,          			  //!<    return reciprocal lattice x index
        int* y,                       //!<    return reciprocal lattice y index
        const int i)                  //!<    lattice index
        const     
    {
        *y = - m_yBandCutOff + i % (2*m_yBandCutOff+1);
        *x = - m_xBandCutOff + floor((double)i/((2*m_yBandCutOff+1)));
        return;
    }
  
    //!
    //!  Map reciprocal lattice indexes x,y, spaning a reciprocal lattice 
    //!  in both +ve and -ve directions, to a single array index i
    //!
    int NonInteractingOflModel::MapLatticeToArray(
        int x,          //!<    return reciprocal lattice x index
        int y) const    //!<    return reciprocal lattice y index
    {
        return (y + m_yBandCutOff)+(2*m_yBandCutOff+1)*(x + m_xBandCutOff);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the total number of terms in the Bloch wave function - this 
    //! will provide the dimension of the difference equation matrix. The number
    //! is given by counting the reciprocal lattice points in a 2D grid 
    //! up to the cut-off size. There is a factor of two to account for the
    //! two spin species. 
    //////////////////////////////////////////////////////////////////////////////////
    iSize_t NonInteractingOflModel::GetDifferenceMatrixDimension() const
    {
        return 2*(2*m_xBandCutOff+1)*(2*m_yBandCutOff+1);
    }

    //!
    //! Get Recoil Energy
    //!
    double NonInteractingOflModel::GetRecoilEnergy() const
    {
        return m_params.m_kappa*m_params.m_kappa/(2.0*m_params.m_mass);
    }

    //!
    //! Return the value of the m_xBandCutOff parameter
    //!
    iSize_t NonInteractingOflModel::GetCutOffX() const
    {
        return m_xBandCutOff;
    }

    //!
    //! Return the value of the m_yBandCutOff parameter
    //!
    iSize_t NonInteractingOflModel::GetCutOffY() const
    {
        return m_yBandCutOff;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Return the nth eigenvector of the difference Hamiltonian, which 
    //! describes the coefficients appearing in the Bloch wave function labelling
    //! each reciprocal lattice point.
    //!
    //! Makes a copy by value of the full eigenvector for the nth band and spin s.
    //! Note that the eigenvector is stored as a staggered array, altenating between 
    //! spin-up and spin down components
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModel::GetBlochCoefficients(
         dcmplx* eigenvector,               //!<    Block of memory to copy the eigenvector into
         const iSize_t n,                   //!<    Band index
         const spin_t s,                    //!<    Spin index (0 or 1) of enum type spin
         const takeConjugate_t conjugate)   //!<    Option to return the complex conjugate or not
         const
    {
        if(n>=0 && n<m_nbrBands)
        {
            if(m_bandsCalculated)
            {
                dcmplx* p_vector = eigenvector;
                dcmplx* p_bloch = m_blochCoefficients;
                p_bloch += n*m_dim + s; 
                for(iSize_t i=0; i<m_dim/2; ++i, ++p_vector, p_bloch+=2)
                {
                    if(_CONJUGATE_ == conjugate)
                    {
                        *p_vector = std::conj(*p_bloch);
                    }
                    else
                    {
                        *p_vector = *p_bloch;
                    }
                }
            }
        }
        else
        {  
            return;
        }    
        return;
    }
    
    //!
    //! Multiply the Bloch vector by a phase factor
    //!
    void NonInteractingOflModel::MultiplyByPhase(const double phaseAngle)
    {
        dcmplx* p_bloch = m_blochCoefficients;
        for(iSize_t i=0; i<m_dim*m_nbrBands; ++i, ++p_bloch)
        {
            *(p_bloch) *= exp(dcmplx(0.0, 1.0)*phaseAngle);
        }
        return;
    }

    //!
    //! Calculate the value of the sigma^z operator applied to the Block
    //! state. This operator will have the effect of muliplying all 
    //! spin-down amplitudes by a minus sign. The spin magnitude is taken to be 1
    //!
    double NonInteractingOflModel::ApplySigmaZ(
        const iSize_t n)           //!<    Band index
        const
    {
        dcmplx* p_bloch = m_blochCoefficients + n*m_dim;
        dcmplx sigmaZ = 0.0;
        for(iSize_t i=0; i<m_dim; ++i, ++p_bloch)
        {
            if(i&1) //  odd terms are associated with spin down
            {
                sigmaZ -= std::conj(*(p_bloch)) * *(p_bloch);
            }
            else
            {
                sigmaZ += std::conj(*(p_bloch)) * *(p_bloch);
            }
        }
        return real(sigmaZ);
    }

    //!
    //! Obtain the value of sum_G a^{kx,ky}_{G,sigma} exp^(ir.(G+{kx,ky})
    //!
    dcmplx NonInteractingOflModel::GetSpatialWaveFunction(
        const spin_t spin,          //!<    Value of the spin index
        const double rx,            //!<    x real-space co-ordinate
        const double ry,            //!<    y real-space co-ordinate
        const double gridSpacing)   //!<    Spacing of the rx,ry grid
        const
    {
        dcmplx* p_bloch  = m_blochCoefficients + spin;
        dcmplx sum = 0.0;
        for(iSize_t i=0; i<m_dim/2; ++i, p_bloch+=2)
        {
            int lCol,mCol;
            this->MapArrayToLattice(&lCol, &mCol, i+spin);
            double kx = m_G.GetKx() + lCol*reciprocalLattice::G1.GetKx() + mCol*reciprocalLattice::G2.GetKx();
            double ky = m_G.GetKy() + lCol*reciprocalLattice::G1.GetKy() + mCol*reciprocalLattice::G2.GetKy();
            sum += *p_bloch * exp(I*(kx*rx*gridSpacing+ky*ry*gridSpacing));
        }
        return sum / sqrt(m_dim);
    }
}   //  End namespace diagonalization 
