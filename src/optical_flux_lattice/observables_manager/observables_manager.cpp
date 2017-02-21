////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines some observables that can be measured for the 
//!     optical flux lattice model 
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
#include "observables_manager.hpp"

namespace diagonalization
{
    //!
    //! Generate a default list of observables in the constructor. 
    //! All observables flags are set to be false by default
    //!
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
        m_observablesList["calculate-translational-density-density"]  = false;
        m_observablesList["calculate-rotational-density-density"] = false; 
        return;
    }

    //!
    //! Copy constructor
    //!
    ObservablesManager::ObservablesManager(ObservablesManager& other)
        :
        m_observablesList(other.m_observablesList),
        m_nbrMostProbable(other.m_nbrMostProbable),
        m_nbrA(other.m_nbrA),
        m_orbitalCut(other.m_orbitalCut),
        m_calcualtedObservables(other.m_calcualtedObservables)
    {
    }

    //!
    //! MPI function to synchronise all variables with those on the given node
    //! 
    void ObservablesManager::MpiSynchronize(
        const int nodeId,                   //!<    Node id to sync with
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        mpi.Sync(&m_nbrMostProbable, 1, nodeId);
        mpi.Sync(&m_nbrA, 1, nodeId);
        mpi.Sync(&m_orbitalCut, 1, nodeId);
        for(const auto &pair : m_observablesList)
        {
            bool temp;
            if(nodeId == mpi.m_id)    //  For the master node
            {
                temp = pair.second;
            }
            mpi.Sync(&temp, 1, nodeId);
            m_observablesList[pair.first] = temp;
        }
    }

    //!
    //! Print a summary of the set observables
    //! 
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
        }
        
        return;
    }
    
    //!
    //! Return the number of observables set
    //!
    int ObservablesManager::GetNbrSetObservables() const
    {
        int counter = 0;
        for(const auto &pair : m_observablesList)
        {
            if(pair.second) ++counter;
        }
        return counter;
    }
    
    //!
    //! Append a observables to the list of observables that are to be 
    //! calculated. Observables are specified by command line arguments
    //!
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
        }
        mpi.ExitFlagTest();
        this->MpiSynchronize(0, mpi);
        return;
    }

    //!
    //! Calculate the values of all the observables in the list
    //! of added observables. If a calculation was unsuccessful, the observables
    //! flag is overwritten as false.
    //!
    void ObservablesManager::CalculateAllObservables(
        OpticalFluxLatticeHamiltonian* hamiltonian, //!<    Address of Hamiltonian class
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
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

    //!
    //! Optionally update the associated SQL database with information on
    //! the calculated observables
    //!
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
                sql.CopyTableToTextFile(sqlTableName,fileName.str());
            }
        } 
        return;
    }
}   //  End namespace diagonalization