////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 25/11/2014
//!                                                                             
//!	 \file
//!     This file defines some observables that can be measured for the 
//!     optical flux lattice model 
//!                                                                  
//!                    Copyright (C) 2014 Simon C Davenport
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

#include "observables_manager.hpp"

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate a default list of observables in the constructor. 
//! All observables flags are set to be false by default
//!
////////////////////////////////////////////////////////////////////////////////

ObservablesManager::ObservablesManager()
    :
    m_nbrMostProbable(0),
    m_nbrA(0),
    m_orbitalCut(0),
    m_calcualtedObservables(false)
{
    m_observablesList["calculate-occupations"]      = false;
    m_observablesList["calculate-susceptibility"]   = false;
    m_observablesList["plot-hamiltonian"]           = false;
    m_observablesList["plot-occupations"]           = false;
    m_observablesList["most-probable-list"]         = false;
    m_observablesList["calculate-density-density"]  = false;
    m_observablesList["calculate-participation-ratio"] = false;
    //m_observablesList["calculate-oes"]              = false;
    //m_observablesList["calculate-pes"]              = false;
    m_observablesList["calculate-translational-density-density"]  = false;
    m_observablesList["calculate-rotational-density-density"] = false;
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Copy constructor
//!
////////////////////////////////////////////////////////////////////////////////

ObservablesManager::ObservablesManager(ObservablesManager& other)
    :
    m_observablesList(other.m_observablesList),
    m_nbrMostProbable(other.m_nbrMostProbable),
    m_nbrA(other.m_nbrA),
    m_orbitalCut(other.m_orbitalCut),
    m_calcualtedObservables(other.m_calcualtedObservables)
{
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief  MPI function to synchronise all variables with those on the given node
//! 
////////////////////////////////////////////////////////////////////////////////

void ObservablesManager::MpiSynchronize(
    const int nodeId,                   //!<    Node id to sync with
    const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{
    mpi.Sync(&m_nbrMostProbable,1,nodeId);
    mpi.Sync(&m_nbrA,1,nodeId);
    mpi.Sync(&m_orbitalCut,1,nodeId);
    
    //  Synchronize the map values over all nodes
    
    for(const auto &pair : m_observablesList)
    {
        bool temp;
        
        if(nodeId == mpi.m_id)    //  For the master node
        {
            temp = pair.second;
        }
    
        mpi.Sync(&temp,1,nodeId);
    
        m_observablesList[pair.first] = temp;
    }
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief  Print a summary of the set observables
//! 
////////////////////////////////////////////////////////////////////////////////

void ObservablesManager::Print(
    const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{
    if(0 == mpi.m_id)    //  For the master node
    {
        utilities::cout.SecondaryOutput()<<"\n\tOBSERVABLES:"<<std::endl;
    
        for(const auto &pair : m_observablesList)
        {
            utilities::cout.SecondaryOutput()<<"\t\t"<<pair.first<<" = "<<pair.second<<std::endl;
        }
        
        if(m_observablesList["most-probable-list"])
        {
            utilities::cout.SecondaryOutput()<<"\t\tNumber of most probable eigs to store = "<<m_nbrMostProbable<<std::endl;
        }

        #if 0
        if(m_observablesList["calculate-pes"])
        {
            utilities::cout.SecondaryOutput()<<"\t\tParticle cut (nA) = "<<m_nbrA<<std::endl;
        }
        
        if(m_observablesList["calculate-oes"])
        {
            utilities::cout.SecondaryOutput()<<"\t\tOrbital cut = "<<m_orbitalCut<<std::endl;
        }
        #endif
    }
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Return the number of observables set
//!
////////////////////////////////////////////////////////////////////////////////

int ObservablesManager::GetNbrSetObservables() const
{
    int counter = 0;

    for(const auto &pair : m_observablesList)
    {
        if(pair.second) counter++;
    }

    return counter;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Append a observables to the list of observables that are to be 
//! calculated. Observables are specified by command line arguments
//!
////////////////////////////////////////////////////////////////////////////////

void ObservablesManager::AddObservables(
    boost::program_options::variables_map* optionList,  //!<    Command line arguments
    utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{
    if(0 == mpi.m_id)    //  For the master node
    {
        for(const auto &pair : m_observablesList)
        {
            m_observablesList[pair.first] = (*optionList)[pair.first].as<bool>();
        }
        
        m_nbrMostProbable = (*optionList)["nbr-most-probable"].as<iSize_t>();
        
        #if 0
        
        m_nbrA = (*optionList)["particle-cut"].as<iSize_t>();
        
        m_orbitalCut = (*optionList)["orbital-cut"].as<iSize_t>();
        
        //  Set default values of orbital and particle cut and check for errors
        
        if(m_nbrA==0)
        {
            //  Set to a half cut
            
            m_nbrA = (*optionList)["nbr"].as<iSize_t>()/2;
        }
        else if(m_nbrA>=(*optionList)["nbr"].as<iSize_t>())
        {
            std::cerr<<"\n\tERROR with particle-cut: must set a value less than the total number of particles"<<std::endl;
            mpi.m_exitFlag = true;
        }
        
        if(m_orbitalCut==0)
        {
            //  Set to a half cut
            
            m_nbrA = ((*optionList)["kx"].as<iSize_t>()*(*optionList)["ky"].as<iSize_t>())/2;
        }
        else if(m_orbitalCut>=((*optionList)["kx"].as<iSize_t>()*(*optionList)["ky"].as<iSize_t>()))
        {
            std::cerr<<"\n\tERROR with orbital-cut: must set a value less than the total number of orbitals"<<std::endl;
            mpi.m_exitFlag = true;
        }
        
        #endif
    }
    
    mpi.ExitFlagTest();

    //  Set variables as on the master node

    this->MpiSynchronize(0,mpi);
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Calculate the values of all the observables in the list
//! of added observables. If a calculation was unsuccessful, the observables
//! flag is overwritten as false.
//!
////////////////////////////////////////////////////////////////////////////////

void ObservablesManager::CalculateAllObservables(
    OpticalFluxLatticeHamiltonian* hamiltonian, //!<    Address of Hamiltonian class
    utilities::MpiWrapper& mpi)           //!<    Instance of the mpi wrapper class
{
    if(m_observablesList["calculate-occupations"])
    {
        m_observablesList["calculate-occupations"] = hamiltonian->CalculateOccupations(mpi);
    }
    
    if(m_observablesList["calculate-susceptibility"])
    {
        m_observablesList["calculate-susceptibility"] = hamiltonian->CalculateSusceptibility(mpi);
    }
    
    if(m_observablesList["plot-hamiltonian"])
    {
        m_observablesList["plot-hamiltonian"] = hamiltonian->PlotHamiltonian(mpi);
    }
    
    if(m_observablesList["most-probable-list"])
    {
        m_observablesList["most-probable-list"] = hamiltonian->ListMostProbable(m_nbrMostProbable,mpi);
    }
    
    if(m_observablesList["calculate-density-density"])
    {
        m_observablesList["calculate-density-density"] = hamiltonian->CalculateDensityDensityFunction(mpi);
    }
    
    if(m_observablesList["calculate-participation-ratio"])
    {
        m_observablesList["calculate-participation-ratio"] = hamiltonian->CalculateParticipationRatio(mpi);
    }
    #if 0
    if(m_observablesList["calculate-oes"])
    {
        m_observablesList["calculate-oes"] = hamiltonian->CalculateOrbitalEntanglementSpectrum(m_orbitalCut,mpi);
    }
    
    if(m_observablesList["calculate-pes"])
    {
        m_observablesList["calculate-pes"] = hamiltonian->CalculateParticleEntanglementSpectrum(m_nbrA,mpi);
    }
    #endif
    if(m_observablesList["calculate-translational-density-density"])
    {
        m_observablesList["calculate-translational-density-density"] = hamiltonian->CalculateTranslationalDensityDensityFunction(mpi);
    }
    
    if(m_observablesList["calculate-rotational-density-density"])
    {
        m_observablesList["calculate-rotational-density-density"] = hamiltonian->CalculateRotationalDensityDensityFunction(mpi);
    }
    
    m_calcualtedObservables = true;
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Optionally update the associated SQL database with information on
//! the calculated observables
//!
////////////////////////////////////////////////////////////////////////////////

void ObservablesManager::UpdateSqlFlags(
    boost::program_options::variables_map* optionList,  //!<    Command line variables map
    const utilities::MpiWrapper& mpi)                   //!<    Instance of the mpi wrapper class
{
    if(0 == mpi.m_id)    //  FOR THE MASTER NODE
    {
        std::string inPath          =  (*optionList)["in-path"].as<std::string>();
        bool useSql                 =  (*optionList)["use-sql"].as<bool>();
        iSize_t sqlId               =  (*optionList)["sql-id"].as<iSize_t>();
        std::string sqlName         =  (*optionList)["sql-name"].as<std::string>();
        std::string sqlTableName    =  (*optionList)["sql-table-name"].as<std::string>();

        if(m_calcualtedObservables && useSql)
        {
            std::stringstream fileName;
                
            fileName.str("");

            fileName<<inPath<<sqlName;

            utilities::Sqlite sql(fileName.str(),_EDIT_EXISTING_);

            //  If any data analysis has been performed at any point then update
            //  the SQL table to indicate that this is so
            
            if(m_observablesList["calculate-occupations"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotOccupations",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }
            
            if(m_observablesList["calculate-susceptibility"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotSusceptibility",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }
            
            if(m_observablesList["most-probable-list"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotMostProbable",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }

            if(m_observablesList["calculate-density-density"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotDensityDensity",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }
            
            if(m_observablesList["calculate-participation-ratio"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotParticipationRatio",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }
            #if 0
            if(m_observablesList["calculate-oes"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotOES",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }
            
            if(m_observablesList["calculate-pes"])
            {
                utilities::SqliteVariable var;
                var.SetValues("GotPES",(int)1);
                sql.UpdateTableEntry(sqlTableName,sqlId,var);
            }
            #endif
            
            //  Update the text copy of the table
                
            sql.CopyTableToTextFile(sqlTableName,fileName.str());
        }
    }
        
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace diagonalization

