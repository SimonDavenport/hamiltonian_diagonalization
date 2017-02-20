////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 13/02/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store a FQHE Haldane pseduopotential
//!     Hamiltonian for in the sphere geometry. 
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

#ifndef _PSEUDOPOTENTIAL_HAMILTONIAN_BASE_HPP_INCLUDED_
#define _PSEUDOPOTENTIAL_HAMILTONIAN_BASE_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

//  Include declaration of associated program options
#include "../program_options/pseudopotential_hamiltonian_options.hpp"

//  Declaration of mState_t, iSize_t and fock_t
#include "../../utilities/general/orbital_and_state_defs.hpp"

//  Data structure containing pseudopotential Hamiltonian data
#include "pseudopotential_hamiltonian_data.hpp"

//  Include various utility classes and definitions
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#include "../../utilities/wrappers/clebsch_gordan_wrapper.hpp"
#include "../../utilities/general/load_bar.hpp"
#include "../../utilities/mathematics/binary_number_tools.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../../utilities/data_structures/sparse_matrix.hpp"
#include "../../hamiltonians/spinless_fermion_hamiltonian.hpp"

#include "angular_momentum_constraint.hpp"

//  Matrix-vector routines
#include "../../hamiltonians/matrix_vector_routines.hpp"

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#include <bitset>
#endif

namespace diagonalization
{

//////      STRUCT/CLASS DECLARATIONS       ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//! \brief The SpherePseudopotentialHamiltonian class defines a FQHE
//! Haldane pseudopotential Hamiltonian the sphere geometry
//!
//////////////////////////////////////////////////////////////////////////////////

template <class H>
class SpherePseudopotentialHamiltonianBase
{
    protected:

    //////      Protected data structures     //////////////////////////////////

    H m_hamiltonian;                    //!<  An object to contain the Hamiltonian itself

    PseudopotentialHamiltonianData m_params;
                                        //!<  Class internal implementation of the     
                                        //!   PseudopotentialHamiltonianData struct      

    //////      Protected functions     ////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Convert signed angular momentum values m into unsigned state
    //! labels k
    //!
    ////////////////////////////////////////////////////////////////////////////////

    kState_t ConvertStateLabel(
        const mState_t m,       //!<    Angular momentum index
        const mState_t maxLz)   //!<   Maximum lz value
        const
    {
        return (-m+maxLz)/2;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Determine the value m1 such that m3+m4-m1-m2 = 0
    //!
    //! \return an integer index defining the m1 orbital, OR an error code
    //! equal to m_maxLz2+2 if momentum conservation is impossible
    //!
    ////////////////////////////////////////////////////////////////////////////////

    mState_t FindM1(
        const mState_t m2,        //!<   Index of m2 vector in the BZ
        const mState_t m3,        //!<   Index of m3 vector in the BZ
        const mState_t m4,        //!<   Index of m4 vector in the BZ
        const mState_t maxLz)     //!<   Maximum lz value
        const
    {
        if(m3+m4-m2 <= maxLz && m3+m4-m2 >= -maxLz)
        {
            return m3+m4-m2;
        }
        else
        {
            //  Return an error code value that's beyond the lz range
            return m_params.m_maxLz2+2;
        }
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate table of 2-body angular momentum conserving sets
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void BaseGenerate4KTable(
        std::vector<kState_t>* kTable,      //!<    Preallocated table of k-values
        const iSize_t llIndex,              //!<   Landau level index
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING MOMENTUM CONSERVING k1,k2,k3,k4 TABLE for LL "<<llIndex<<std::endl;
        }

        mState_t maxLz;

        if(llIndex == 0)
        {
            maxLz = m_params.m_maxLz;
        }
        else if(llIndex == 1)
        {
            maxLz = m_params.m_maxLz2;
        }
        else
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                utilities::cout.AdditionalInfo()<<"\n\tERROR in BaseGenerate4KTable: Invalid ll index "<<llIndex<<std::endl;
            }
            
            exit(EXIT_FAILURE);
        }

        auto it_table = kTable->begin();

        for(mState_t m4=maxLz;m4>=-maxLz;m4-=2)
        {
            for(mState_t m3=maxLz;m3>=-maxLz;m3-=2)
            {
                for(mState_t m2=maxLz;m2>=-maxLz;m2-=2,++it_table)
                {
                    if(it_table < kTable->end())
                    {
                        *(it_table) = this->ConvertStateLabel(this->FindM1(m2,m3,m4,maxLz),maxLz);
                    }
                }
            }
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Convert two-body pseudopotential coefficients into a list of
    //! matrix elements Vkkkk.
    //!
    //! The coefficients are given by:
    //!
    //! \[ V_{k1,k2,k3,k4}  = \delta (k3+k4-k1-k2) * sum_{L=0^m_maxLz} sum_{M=-L^L}
    //!                                             V_L <L,M||k1>|k2> <k4|<k3||L,M>
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void BaseGetCoefficientsFromPseudopotentials(
        std::vector<double>* listOfQuarticTerms,
                                            //!<   List of coefficients to populate
        const iSize_t llIndex,              //!<   Landau level index
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING HAMILTONIAN COEFFICIENTS FROM TWO-BODY PSEUDOPOTENTIALS for LL "<<llIndex<<std::endl;
        }
    
        mState_t maxLz;
        std::vector<double> pseudopotentials;
    
        if(llIndex == 0)
        {
            maxLz = m_params.m_maxLz;
            pseudopotentials = m_params.m_twoBodyPseudopotentials;
        }
        else if(llIndex == 1)
        {
            maxLz = m_params.m_maxLz2;
            pseudopotentials = m_params.m_twoBodyPseudopotentials2LL;
        }
        else
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                utilities::cout.AdditionalInfo()<<"\n\tERROR in GetCoefficientsFromPseudopotentials: Invalid ll index "<<llIndex<<std::endl;
            }
            
            exit(EXIT_FAILURE);
        }
    
        auto it_table = listOfQuarticTerms->begin();

        //  Generate a look-up table of Clebsch Gordan coefficients
        
        utilities::ClebschGordanCoefficients<double> clebschGordan;

        clebschGordan.AddMultiplet(maxLz,maxLz);

        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATED TEMPORARY TABLE OF "<<clebschGordan.GetSize()<<" CLEBSCH GORDAN COEFFICIENTS."<<std::endl;
        }
        
        for(mState_t m4=maxLz;m4>=-maxLz;m4-=2)
        {
            for(mState_t m3=maxLz;m3>=-maxLz;m3-=2)
            {
                for(mState_t m2=maxLz;m2>=-maxLz;m2-=2,++it_table)
                {
                    if(it_table != listOfQuarticTerms->end())
                    {
                        double total = 0.0;

                        mState_t m1 = this->FindM1(m2,m3,m4,maxLz);
                        
                        if(m1 <= maxLz)
                        {
                            //  Perform sum over L and M
                            
                            for(mState_t L=0;L<2*pseudopotentials.size();L+=2)
                            {
                                double temp = 0.0;
                                
                                //  For consistency with known results, the L value must
                                //  be reflected in this way:
                                
                                mState_t Lact = 2*maxLz-L;
                                mState_t M = m1+m2;
                                
                                //std::cout<<"\t\tL = "<<(int)Lact<<std::endl;
                                
                                if(-Lact<=M && M<=Lact)
                                {
                                    //  Add terms of the form
                                    //  <L,M||m1>|m2> <m4|<m3||L,M>

                                    //temp += random();
                                    
                                    temp += clebschGordan.Value(maxLz,m1,maxLz,m2,Lact,M)*clebschGordan.Value(maxLz,m3,maxLz,m4,Lact,M);
                                    
                                }
                                
                                if(L/2<pseudopotentials.size())
                                {
                                    total -= pseudopotentials[L/2]*temp;
                                }
                            }
                            
                            #if _DEBUG_
                            if(fabs(total)>0.0000000001)    
                            {
                                std::cout<<"\t<"<<m1<<" "<<m2<<"|V|"<<m3<<" "<<m4<<"> = "<<total<<std::endl;
                            }
                            #endif
                        }
                        
                        *it_table = total;
                    }
                }
            }
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Diagonalize our Hamiltonian using a selected method
    //! 
    ////////////////////////////////////////////////////////////////////////////////

    template <class L1,class L2>
    void BaseDiagonalize(
        L1* quadraticLookUpTables,          //!<    Class container for look-up tables
        L2* quarticLookUpTables,            //!<    Class container for look-up tables
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if(m_params.m_hamiltonianBuilt)
        {
            //////  CALL SELECTED DIAGONALIZATION ROUTINE

            //  Switch to LAPACK for very small matrices

            if(m_hamiltonian.m_data.m_fockSpaceDim<60)
            {
                if(0 == mpi.m_id)
                {
                    utilities::cout.SecondaryOutput()<<"\n\tMATRIX DIMENSION <60 - SETTING DIAGONALIZATION METHOD TO FULL\n"<<std::endl;
                }
                
                m_params.m_method = _FULL_;
            }

            if(_FULL_ == m_params.m_method)
            {
                m_params.m_hamiltonianDiagonalized = m_hamiltonian.FullDiagonalize(m_params.m_nbrEigenvaluesToFind,mpi);
            }
            else if(_LANCZOS_ == m_params.m_method)
            {
                //  Update the file names for the input/output ARPACK vectors
                
                std::stringstream inFileName;
                std::stringstream outFileName;
                
                inFileName.str("");
                outFileName.str("");
                
                inFileName<<m_params.m_inPath<<"/"<<m_params.m_initialVectorFile;
                outFileName<<m_params.m_outPath<<"/"<<m_params.m_finalVectorFile;
                
                if(m_params.m_blockDiagonalize)
                {
                    inFileName<<"_sector_"<<m_params.m_totalLz;
                    outFileName<<"_sector_"<<m_params.m_totalLz;
                }
                
                inFileName<<"_proc_"<<mpi.m_id<<".bin";
                outFileName<<"_proc_"<<mpi.m_id<<".bin";
                
                m_hamiltonian.UpdateArpackFileNames(inFileName.str(),outFileName.str(),mpi);
            
                MatrixVectorFunction<double,L1,L2> mv(m_hamiltonian,mpi);

                m_params.m_hamiltonianDiagonalized = m_hamiltonian.LanczosDiagonalize(
                quadraticLookUpTables,quarticLookUpTables,mv,
                m_params.m_nbrEigenvaluesToFind,mpi);
            }

            if(!m_params.m_hamiltonianDiagonalized)
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::cerr<<"\n\tWARNING - NOT ABLE TO RUN DIAGONALIZATION ROUTINE."<<std::endl;
                }
            }
        }
    }

    //////      Public interface functions     ////////////////////////////////

    public:
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Default constructor for the SpherePseudopotentialHamiltonian
    //! class
    //!
    ////////////////////////////////////////////////////////////////////////////////

    SpherePseudopotentialHamiltonianBase()
    {
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor from command line arguments
    //!
    ////////////////////////////////////////////////////////////////////////////////

    SpherePseudopotentialHamiltonianBase(
        boost::program_options::variables_map* optionList,    
                                    //!<    Parsed command line argument list
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    {
        m_params = PseudopotentialHamiltonianData(optionList,mpi);

        //  Set ARPACK related options
        
        m_hamiltonian.InitializeArpack(optionList,mpi);
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Destructor
    //! 
    ////////////////////////////////////////////////////////////////////////////////

    ~SpherePseudopotentialHamiltonianBase()
    {
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Set total lz sector
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void SetSector(
        const mState_t lzSector)     //!<    Specify total angular momentum sector
    {
        m_params.m_totalLz = lzSector;
        m_params.m_blockDiagonalize = true;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Set the number of eigenvalues to find
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void SetNbrEigenvaluesToFind(
        const iSize_t nbrEigenvaluesToFind)     //!<    Specify new number of eigenvalues to find
    {
        m_params.m_nbrEigenvaluesToFind = nbrEigenvaluesToFind;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////      Pure virtual functions      ////////////////////////////////////////////////////////
    
    //  Build look-up tables
    virtual void BuildLookUpTables(const utilities::MpiWrapper& mpi)=0;

    // Build Fock basis
    virtual void BuildFockBasis(utilities::MpiWrapper& mpi)=0;
    
    //  Put entries into full Hamiltonian
    virtual void BuildHamiltonian(utilities::MpiWrapper& mpi)=0;
    
    //  Diagonalize the Hamiltonian 
    virtual void Diagonalize(utilities::MpiWrapper& mpi)=0;
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to clear any currently stored FermionHamiltonian object and
    //! reset the m_params.m_hamiltonianDiagonalized and m_params.m_hamiltonianBuilt flags
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void ClearHamiltonian()
    {
        if(m_params.m_hamiltonianBuilt)
        {
            m_hamiltonian.Clear();
            m_params.m_hamiltonianBuilt = false;
            m_params.m_hamiltonianDiagonalized = false;
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Make a copy of the eigenvalue data to an external array
    //!
    ////////////////////////////////////////////////////////////////////////////////
  
    void GetEigenvalues(
        double* buffer,         //!<    Buffer to store returned eigenvalues
        const iSize_t nbrEigenvalues)
                                //!<    Number of eigenvalues
        const
    {
        m_hamiltonian.GetEigenvalues(buffer,nbrEigenvalues);
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Output all stored eigenvalue and eigenvector data to a text file
    //!
    //! This function MUST be called on all nodes if run in parallel 
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void EigensystemToFile(
        const bool writeEigenvalues,    //!<    Option to write eigenvalues to file
        const bool writeEigenvectors,   //!<    Option to write eigenvectors to file
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {   
        if(writeEigenvalues || (writeEigenvectors && m_params.m_hamiltonianDiagonalized))
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t============ WRITING EIGENSYSTEM DATA TO FILE ============ "<<std::endl;
            }
        
            std::stringstream ss;
            
            ss.str("");
            
            ss<<m_params.MakeBaseFileName(_OUT_)<<"_eigensystem.dat";
            
            std::vector<iSize_t> permsList(m_params.m_nbrEigenvaluesToFind);
            
            for(iSize_t i=0;i<m_params.m_nbrEigenvaluesToFind;++i)
            {
                permsList[i] = i;
            }
            
            m_hamiltonian.EigensystemToFile(ss.str(),writeEigenvalues,writeEigenvectors,mpi);

            mpi.ExitFlagTest();
        }

        return;
    }   

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Obtain eigenvalue and eigenvector data from existing files
    //!
    //! This function MUST be called on all nodes if run in parallel 
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void EigensystemFromFile(
        const bool readEigenvalues,
        const bool readEigenvectors,
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {   
        if((readEigenvalues || readEigenvectors) && !m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {   
                utilities::cout.SecondaryOutput()<<"\n\t============ READING EIGENSYSTEM DATA FROM FILE ============ "<<std::endl;
            }
        
            m_hamiltonian.Initialize(m_params.m_nbrParticles,m_params.m_nbrOrbitals,utilities::_SPARSE_MAPPED_,mpi);
        
            //  If the matrix is Block-diagonalized then generate an appropriate set
            //  of allowed Fock basis states
        
            if(m_params.m_blockDiagonalize)
            {
                //  For an explanation of this function call see the comments
                //  in the BuildHamiltonian() function
                
                FermionFockBasis basis;
                
                basis.GenerateFockSpace(m_params.m_nbrParticles,m_params.m_nbrOrbitals,std::bind(&TestAngularMomentumSector,std::placeholders::_1,m_params.m_totalLz,m_params.m_nbrOrbitals),mpi);
            
                m_hamiltonian.SetFockBasis(basis);
            }
        
            std::stringstream ss;
            
            ss.str("");
            
            ss<<m_params.MakeBaseFileName(_IN_)<<"_eigensystem.dat";

            m_hamiltonian.EigensystemFromFile(ss.str(),readEigenvalues,readEigenvectors,mpi);
            
            mpi.ExitFlagTest();
            
            //  Counts as if the Hamiltonian were diagonalized

            m_params.m_hamiltonianDiagonalized = true;
        }

        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Write currently stored Hamiltonian matrix data to a file
    //!
    ////////////////////////////////////////////////////////////////////////////////

    void HamiltonianToFile(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(m_params.m_hamiltonianBuilt)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t============ WRITING HAMILTONIAN DATA TO FILE ============ "<<std::endl;
            }

            std::stringstream ss;

            ss.str("");

            ss<<m_params.MakeBaseFileName(_OUT_)<<"_hamiltonian_"<<mpi.m_id<<".dat";

            m_hamiltonian.HamiltonianToFile(ss.str(),mpi);

            mpi.ExitFlagTest();
        }
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

};

}   //  End namespace diagonalization

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#endif

