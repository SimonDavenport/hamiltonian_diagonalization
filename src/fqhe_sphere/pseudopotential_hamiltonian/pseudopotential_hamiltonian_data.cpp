////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 13/02/2015
//!                                                                             
//!	 \file
//!     This file defines a class to store a parameters defining a
//!     FQHE Haldane pseduopotential Hamiltonian in the sphere geometry
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
#include "pseudopotential_hamiltonian_data.hpp"

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Default Constructor
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData::PseudopotentialHamiltonianData()
:   m_nbrOrbitals(9),
    m_nbrParticles(3),
    m_maxLz(8),
    m_maxLz2(10),
    m_totalLz(0),
    m_nbrLevels(1),
    m_outPath("./"),
    m_inPath("./"),
    m_outFileName("out"),
    m_nbrEigenvaluesToFind(4),
	m_blockDiagonalize(false),
    m_method(_FULL_),
    m_initialVectorFile("initVector.bin"),
    m_finalVectorFile("finalVector.bin"),
    m_lookupTablesBuilt(false),
    m_fockBasisBuilt(false),
    m_hamiltonianBuilt(false),
    m_hamiltonianDiagonalized(false)
{
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Constructor from minimal data for single LLL problem
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData::PseudopotentialHamiltonianData(
    const iSize_t nbrParticles,                     //!<    Number of particles
    const iSize_t nbrOrbitals,                      //!<    Number of orbitals
    const std::vector<double>& pseudopotentials)    //!<    List of two-body pseudopotentials
:   m_nbrOrbitals(nbrOrbitals),
    m_nbrParticles(nbrParticles),
    m_totalLz(0),
    m_nbrLevels(1),
    m_twoBodyPseudopotentials(pseudopotentials),
    m_outPath("./"),
    m_inPath("./"),
    m_outFileName("out"),
    m_nbrEigenvaluesToFind(4),
	m_blockDiagonalize(false),
    m_method(_FULL_),
    m_initialVectorFile("initVector.bin"),
    m_finalVectorFile("finalVector.bin"),
    m_lookupTablesBuilt(false),
    m_fockBasisBuilt(false),
    m_hamiltonianBuilt(false),
    m_hamiltonianDiagonalized(false)
{
    this->SetMaxLz(m_nbrOrbitals,m_nbrLevels);
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Constructor from minimal data including 2LLs
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData::PseudopotentialHamiltonianData(
    const iSize_t nbrParticles,                     //!<    Number of particles
    const iSize_t nbrOrbitals,                      //!<    Number of orbitals
    const std::vector<double>& pseudopotentials,    //!<    List of two-body pseudopotentials
    const std::vector<double>& pseudopotentials2LL) //!<    List of 2LL two-body pseudopotentials
:   m_nbrOrbitals(nbrOrbitals),
    m_nbrParticles(nbrParticles),
    m_totalLz(0),
    m_nbrLevels(2),
    m_twoBodyPseudopotentials(pseudopotentials),
    m_twoBodyPseudopotentials2LL(pseudopotentials2LL),
    m_outPath("./"),
    m_inPath("./"),
    m_outFileName("out"),
    m_nbrEigenvaluesToFind(4),
	m_blockDiagonalize(false),
    m_method(_FULL_),
    m_initialVectorFile("initVector.bin"),
    m_finalVectorFile("finalVector.bin"),
    m_lookupTablesBuilt(false),
    m_fockBasisBuilt(false),
    m_hamiltonianBuilt(false),
    m_hamiltonianDiagonalized(false)
{
    this->SetMaxLz(m_nbrOrbitals,m_nbrLevels);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
////////////////////////////////////////////////////////////////////////////////
//! \brief Destructor
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData::~PseudopotentialHamiltonianData()
{
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
////////////////////////////////////////////////////////////////////////////////
//! \brief Copy constructor
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData::PseudopotentialHamiltonianData(
    const PseudopotentialHamiltonianData& other) //!<    Instance of the object
:   m_nbrOrbitals(other.m_nbrOrbitals),
    m_nbrParticles(other.m_nbrParticles),
    m_maxLz(other.m_maxLz),
    m_maxLz2(other.m_maxLz2),
    m_totalLz(other.m_totalLz),
    m_nbrLevels(other.m_nbrLevels),
    m_twoBodyPseudopotentials(other.m_twoBodyPseudopotentials),
    m_twoBodyPseudopotentials2LL(other.m_twoBodyPseudopotentials2LL),
    m_outPath(other.m_outPath),
    m_inPath(other.m_inPath),
    m_outFileName(other.m_outFileName),
    m_nbrEigenvaluesToFind(other.m_nbrEigenvaluesToFind),
	m_blockDiagonalize(other.m_blockDiagonalize),
    m_method(other.m_method),
    m_initialVectorFile(other.m_initialVectorFile),
    m_finalVectorFile(other.m_finalVectorFile),
    m_lookupTablesBuilt(other.m_lookupTablesBuilt),
    m_fockBasisBuilt(other.m_fockBasisBuilt),
    m_hamiltonianBuilt(other.m_hamiltonianBuilt),
    m_hamiltonianDiagonalized(other.m_hamiltonianDiagonalized)
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Overload assignment operator
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData& PseudopotentialHamiltonianData::operator=(
    const PseudopotentialHamiltonianData& other)   //!<    Instance of the object
{
    m_nbrOrbitals = other.m_nbrOrbitals;
    m_nbrParticles = other.m_nbrParticles;
    m_maxLz = other.m_maxLz;
    m_maxLz2 = other.m_maxLz2;
    m_nbrLevels = other.m_nbrLevels;
    m_totalLz = other.m_totalLz;
    m_twoBodyPseudopotentials = other.m_twoBodyPseudopotentials;
    m_twoBodyPseudopotentials2LL = other.m_twoBodyPseudopotentials2LL;
    m_outPath = other.m_outPath;
    m_inPath = other.m_inPath;
    m_outFileName = other.m_outFileName;
    m_nbrEigenvaluesToFind = other.m_nbrEigenvaluesToFind;
	m_blockDiagonalize = other.m_blockDiagonalize;
    m_method = other.m_method;
    m_initialVectorFile = other.m_initialVectorFile;
    m_finalVectorFile = other.m_finalVectorFile;
    m_lookupTablesBuilt = other.m_lookupTablesBuilt;
    m_fockBasisBuilt = other.m_fockBasisBuilt;
    m_hamiltonianBuilt = other.m_hamiltonianBuilt;
    m_hamiltonianDiagonalized = other.m_hamiltonianDiagonalized;

    return *this;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief  MPI function to synchronise all variables with those on the given node
//! 
////////////////////////////////////////////////////////////////////////////////

void PseudopotentialHamiltonianData::MpiSynchronize(
    const int nodeId,          //!<    Node id to sync with
    const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{ 
    mpi.Sync(&m_nbrOrbitals,1,nodeId);
    mpi.Sync(&m_nbrParticles,1,nodeId);
    mpi.Sync(&m_maxLz,1,nodeId);
    mpi.Sync(&m_maxLz2,1,nodeId);
    mpi.Sync(&m_totalLz,1,nodeId);
    mpi.Sync(&m_nbrLevels,1,nodeId);
    mpi.Sync(&m_twoBodyPseudopotentials,nodeId);
    mpi.Sync(&m_twoBodyPseudopotentials2LL,nodeId);
    mpi.Sync(m_inPath,nodeId);
    mpi.Sync(m_outPath,nodeId);
    mpi.Sync(m_outFileName,nodeId);
    mpi.Sync(&m_nbrEigenvaluesToFind,1,nodeId);
    mpi.Sync(&m_blockDiagonalize,1,nodeId);
    
    int method = 0 ;
    if(_FULL_ == m_method)
    {
        method = 0;
    }
    else if(_LANCZOS_ == m_method)
    {
        method = 1;
    }

    mpi.Sync(&method,1,nodeId);
    
    if(0 == method)
    {
        m_method = _FULL_;
    }
    else if(1 == method)
    {
        m_method = _LANCZOS_;
    }

    mpi.Sync(m_initialVectorFile,nodeId);
    mpi.Sync(m_finalVectorFile,nodeId);
    mpi.Sync(&m_lookupTablesBuilt,1,nodeId);
    mpi.Sync(&m_fockBasisBuilt,1,nodeId);
    mpi.Sync(&m_hamiltonianBuilt,1,nodeId);
    mpi.Sync(&m_hamiltonianDiagonalized,1,nodeId);

    //	Barrier waits for all nodes to catch up to this point
	MPI_Barrier(mpi.m_comm);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Constructor for the PseudopotentialHamiltonianData struct
//! from command line arguments
//!
////////////////////////////////////////////////////////////////////////////////

PseudopotentialHamiltonianData::PseudopotentialHamiltonianData(
    boost::program_options::variables_map* optionList,  
                                //!<    Command line arguments list
    utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    :
    m_totalLz(0),
	m_lookupTablesBuilt(false),
    m_hamiltonianBuilt(false),
    m_hamiltonianDiagonalized(false)
{
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    //////      Read in parameter values from command line parser       ////
	
	    m_nbrOrbitals           =  (*optionList)["nbr-orbitals"].as<iSize_t>();
	    m_nbrParticles          =  (*optionList)["nbr-particles"].as<iSize_t>();
        m_outPath               =  (*optionList)["out-path"].as<std::string>();
        m_inPath                =  (*optionList)["in-path"].as<std::string>();
        m_nbrEigenvaluesToFind  =  (*optionList)["nbr-eigenvalues"].as<iSize_t>();
        m_initialVectorFile     =  (*optionList)["initial-file"].as<std::string>();
        m_finalVectorFile       =  (*optionList)["final-file"].as<std::string>();
        m_blockDiagonalize      =  (*optionList)["block-diagonalize"].as<bool>();
        m_nbrLevels             =  (*optionList)["nbr-levels"].as<iSize_t>();

        if(optionList->count("two-body-pseudopotentials"))
        {
            m_twoBodyPseudopotentials = (*optionList)["two-body-pseudopotentials"].as<std::vector<double> >();
        }
        else
        {
            m_twoBodyPseudopotentials.push_back(0);
            m_twoBodyPseudopotentials.push_back(1);
        }
        
        if(optionList->count("two-body-pseudopotentials-2ll"))
        {
            m_twoBodyPseudopotentials2LL = (*optionList)["two-body-pseudopotentials-2ll"].as<std::vector<double> >();
        }
        else
        {
            m_twoBodyPseudopotentials2LL.push_back(0);
            m_twoBodyPseudopotentials2LL.push_back(1);
        } 
                
        //////      Perform consistency checks       ///////////////////////////
        {
            if(m_nbrLevels==2 && (2 >= m_nbrOrbitals || (m_nbrOrbitals & 1)))
            {
                std::cerr<<"\n\tERROR: MUST SPECIFY nbr orbitals > 2 and an even number of them."<<std::endl;
                mpi.m_exitFlag = true;
                goto errorPoint;
            }
            
            if(m_nbrLevels>2)
            {
                std::cerr<<"\n\tERROR: cannot specify more than 2 Landau levels"<<std::endl;
                mpi.m_exitFlag = true;
                goto errorPoint;
            }

	        this->SetMaxLz(m_nbrOrbitals,m_nbrLevels);
	        
	        //////      Check that the number of pseudopotentials       ////////////
            //////      specified does not exceed 2S                    ////////////
            
            if((int)m_twoBodyPseudopotentials.size()>m_maxLz)
            {
                std::cerr<<"\n\tERROR: TOO MANY 2-BODY PSEUDOPOTENTIAL COEFFICIENTS SPECIFIED"<<std::endl;
                mpi.m_exitFlag = true;
                goto errorPoint;
            }
            
            if((int)m_twoBodyPseudopotentials2LL.size()>m_maxLz2)
            {
                std::cerr<<"\n\tERROR: TOO MANY 2-BODY 2LL PSEUDOPOTENTIAL COEFFICIENTS SPECIFIED"<<std::endl;
                mpi.m_exitFlag = true;
                goto errorPoint;
            }

	        //////      Print out parameter values      ////////////////////////////
	        
	        utilities::cout.MainOutput()<<"\n\tPSEUDOPOTENTIAL HAMILTONIAN PARAMETERS:"<<std::endl;
            utilities::cout.MainOutput()<<"\n\t\tNumber of orbitals: \t\t"<<m_nbrOrbitals<<std::endl;
            utilities::cout.MainOutput()<<"\t\tNumber of LLs: \t\t\t"<<m_nbrLevels<<std::endl;
            utilities::cout.MainOutput()<<"\t\tMonopole strength 2S: \t\t"<<m_nbrOrbitals-1<<std::endl;
            utilities::cout.MainOutput()<<"\t\tNumber of particles: \t\t"<<m_nbrParticles<<std::endl;
            utilities::cout.MainOutput()<<"\t\tNumber of eigs to find: \t"<<m_nbrEigenvaluesToFind<<std::endl;
            utilities::cout.MainOutput()<<"\t\tOut Path: \t\t"<<m_outPath<<std::endl;
            utilities::cout.MainOutput()<<"\t\tIn Path: \t\t"<<m_inPath<<std::endl;   
            
            utilities::cout.SecondaryOutput()<<"\n\tORBITAL LABELS:"<<std::endl<<"\t";

            if(m_nbrLevels==1)
	        {
                for(mState_t i=m_maxLz;i>=-m_maxLz;i-=2)
	            {   
	                if(i & 1)
	                {
                        utilities::cout.SecondaryOutput()<<"\t"<<i;
                        
                        if(0!=i)    utilities::cout.SecondaryOutput()<<"/2";
                    }
                    else
                    {
                        utilities::cout.SecondaryOutput()<<"\t"<<i/2;
                    }
                }
            }
            
            if(m_nbrLevels==2)
	        {
	            for(mState_t i=m_maxLz2;i>=-m_maxLz2;i-=2)
	            {   
	                if(i & 1)
	                {
                        utilities::cout.SecondaryOutput()<<"\t"<<i;
                        
                        if(0!=i)    utilities::cout.SecondaryOutput()<<"/2";
                    }
                    else
                    {
                        utilities::cout.SecondaryOutput()<<"\t"<<i/2;
                    }
                }
	        
	            utilities::cout.SecondaryOutput()<<std::endl<<"\t\t";

	            for(mState_t i=m_maxLz;i>=-m_maxLz;i-=2)
	            {   
	                if(i & 1)
	                {
                        utilities::cout.SecondaryOutput()<<"\t"<<i;
                        
                        if(0!=i)    utilities::cout.SecondaryOutput()<<"/2";
                    }
                    else
                    {
                        utilities::cout.SecondaryOutput()<<"\t"<<i/2;
                    }
                }
	        }
            
            utilities::cout.SecondaryOutput()<<std::endl;
            
	        int tmp = (*optionList)["method"].as<iSize_t>();
            
            if(0==tmp)	    m_method = _FULL_;
            else if(1==tmp)	m_method = _LANCZOS_;
            else 
            {
                std::cerr<<"\n\tERROR: INVALID METHOD : "<<tmp<<std::endl;
                mpi.m_exitFlag = true;
            }

            if(_FULL_ == m_method)
            {
                utilities::cout.SecondaryOutput()<<"\n\tDiagonalization method: Full (LAPACK)"<<std::endl;
            }
            else if(_LANCZOS_ == m_method)
            {
                utilities::cout.SecondaryOutput()<<"\n\tDiagonalization method: Lanczos (ARPACK)"<<std::endl;
            }

            utilities::cout.SecondaryOutput()<<"\n\tPSEUDOPOTENTIAL COEFFICIENTS:\n"<<std::endl;

	        for(int i=0;i<(int)m_twoBodyPseudopotentials.size();++i)
	        {
	            utilities::cout.SecondaryOutput()<<"\t\tV_"<<i<<" = "<<m_twoBodyPseudopotentials[i]<<std::endl;
	        }

            utilities::cout.SecondaryOutput()<<"\n\t2ND LL PSEUDOPOTENTIAL COEFFICIENTS: (USED IF ENABLED)\n"<<std::endl;

            utilities::cout.SecondaryOutput()<<std::endl;
        
            for(int i=0;i<(int)m_twoBodyPseudopotentials2LL.size();++i)
            {
                utilities::cout.SecondaryOutput()<<"\t\tV^2_"<<i<<" = "<<m_twoBodyPseudopotentials2LL[i]<<std::endl;
            }
        
        }
        
        errorPoint:

        std::stringstream fileName;
    
        fileName.str("");
        
        fileName<<"fqhe_sphere_n_"<<m_nbrParticles<<"_2s_"<<m_maxLz;
        
        m_outFileName = fileName.str();
	}
	
	mpi.ExitFlagTest();

	//  Sync the parameter values with those set on the master node
	
	this->MpiSynchronize(0,mpi);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate an output file name base, using struct parameters
//!
////////////////////////////////////////////////////////////////////////////////

std::string PseudopotentialHamiltonianData::MakeBaseFileName(
    const io_t io)  //!<    Flag to set input/output file name
const
{
    std::stringstream fileName;
            
    fileName.str("");
    
    if(_OUT_ == io)
    {
        fileName<<m_outPath;
    }
    else if(_IN_ == io)
    {
        fileName<<m_inPath;
    }
    
    fileName<<"/"<<m_outFileName;
    
    if(m_blockDiagonalize)
    {
        fileName<<"_sector_"<<m_totalLz;
    }
    
    return fileName.str();
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Set the minimum value of Lz such that the orbitals are labelled by 
//! -m_maxLz to m_maxLz in steps of 2
//!
////////////////////////////////////////////////////////////////////////////////

void PseudopotentialHamiltonianData::SetMaxLz(
    const iSize_t nbrOrbitals,  //!<    Number of orbitals
    const iSize_t nbrLevels)    //!<    Number of Landau levels
{
    if(nbrLevels==1)
    {
        //  Label each level in the same way and treat independently
        m_maxLz  = m_nbrOrbitals-1;
        m_maxLz2 = m_nbrOrbitals-1;
    }
    else if(nbrLevels==2)
    {   
        //  Line up the angular momentum labels of the two levels
        iSize_t nbrOrbitalsLLL = (m_nbrOrbitals-2)/2;

        m_maxLz = nbrOrbitalsLLL-1;
        m_maxLz2 = nbrOrbitalsLLL+1;
    }
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace diagonalization
