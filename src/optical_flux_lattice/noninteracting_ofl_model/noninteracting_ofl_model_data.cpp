////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                 
//!                                                                             
//!	 \file
//!     This file defines a class to store data describing the non-interacting
//!     optical flux lattice model. See e.g. PRL 109, 265301 (2012).
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
#include "noninteracting_ofl_model_data.hpp"

namespace diagonalization
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Default Constructor
    //!
    //!  Set default model parameters (roughly corresponding to the flattest
    //!  band part of the phase diagram)
    ////////////////////////////////////////////////////////////////////////////////
    NonInteractingOflModelData::NonInteractingOflModelData()
    :   m_theta(0.3),
        m_V0(4.0),
        m_epsilon(0.4),
        m_kappa(1.0),
        m_mass(1.0)
    {}

    //!
    //! Copy Constructor
    //!
    NonInteractingOflModelData::NonInteractingOflModelData(
        const NonInteractingOflModelData& other)  //!<    Instance of the object
    :   m_theta(other.m_theta),
        m_V0(other.m_V0),
        m_epsilon(other.m_epsilon),
        m_kappa(other.m_kappa),
        m_mass(other.m_mass)
    {}

    //!
    //! Assignment operator overload
    //!
    NonInteractingOflModelData& NonInteractingOflModelData::operator=(
        const NonInteractingOflModelData& other)  //!<    Instance of the object
    {
        m_theta = other.m_theta;
        m_V0 = other.m_V0;
        m_epsilon = other.m_epsilon;
        m_kappa = other.m_kappa;
        m_mass = other.m_mass;
        return *this;
    }

    //!
    //! A constructor for the Hamiltonian parameters stuct that gets the
    //! parameter values from a file defined by the command line arguments
    //!
    NonInteractingOflModelData::NonInteractingOflModelData(
        boost::program_options::variables_map* optionList,  
                                    //!<    Parsed command line argument list
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            std::stringstream fileName;
	        fileName.str("");
            bool useSql = (*optionList)["use-sql"].as<bool>();
            if(useSql)
            {
                //  If sql option is set, then look for model data in an sqlite file
                fileName<<(*optionList)["in-path"].as<std::string>()<<(*optionList)["sql-name"].as<std::string>();
                this->ReadFromSql((*optionList)["sql-table-name"].as<std::string>(), fileName.str(), (*optionList)["sql-id"].as<iSize_t>(), mpi);
                if(mpi.m_exitFlag) return;
            }
            else
            { 
                //  Default to look for model data in the specified text file
                fileName<<(*optionList)["in-path"].as<std::string>()<<(*optionList)["params-file"].as<std::string>();
                this->ReadFromFile(fileName.str(),mpi);
                if(mpi.m_exitFlag) return;
            }
            //  Print out summary of model parameters:
            utilities::cout.MainOutput()<<"\n\tSINGLE PARTICLE MODEL PARAMETERS:\n"<<std::endl;
            utilities::cout.MainOutput()<<"\t\ttheta: \t\t"<<m_theta<<std::endl;
            utilities::cout.MainOutput()<<"\t\tV0: \t\t"<<m_V0<<std::endl;
            utilities::cout.MainOutput()<<"\t\tepsilon: \t"<<m_epsilon<<std::endl;
            utilities::cout.MainOutput()<<"\t\tkappa: \t\t"<<m_kappa<<std::endl;
            utilities::cout.MainOutput()<<"\t\tmass: \t\t"<<m_mass<<std::endl;
        }
        mpi.ExitFlagTest();
	    this->MpiSynchronize(0, mpi);
    }

    //!
    //! Read model parameters for a specified file.
    //! 
    void NonInteractingOflModelData::ReadFromFile(
        std::string fileName,           //!<    Name of file where parameters are contained
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        std::ifstream f_param;
        std::string format = "text";
        utilities::GenFileStream(f_param, fileName, format, mpi);
	    if(mpi.m_exitFlag)
	    {
		    return;
	    }
	    else
	    {
	        utilities::cout.MainOutput()<<"\n\tReading model parameters from file "<<fileName<<std::endl;
	    }
	    std::string line;
	    while(getline(f_param,line))
	    {
		    if(line[0]!='#')	break;
	    }
	    f_param >> m_theta >> m_V0 >> m_epsilon >> m_kappa >> m_mass;
        f_param.close();
        return;
    }

    //!
    //! Generate an output file name containing the current model parameters
    //!
    void NonInteractingOflModelData::GenerateParametersFile() const
    {
        std::stringstream parameters;
        parameters<<
        "####################################################\n"
        "##    This file contains a set of single particle   \n"
        "##    Hamiltonian parameters for the optical flux   \n"
        "##    lattice model.                              \n\n"
        <<m_theta<<"\n"
	    <<m_V0<<"\n"
	    <<m_epsilon<<"\n"
	    <<m_kappa<<"\n"
	    <<m_mass<<"\n"; 
        std::ofstream f_parameters;
        f_parameters.open("./parameters.dat", std::ios::out);
        f_parameters<<parameters.str().c_str();
        f_parameters.close();
    }

    //!
    //! Function to read parameters from a sql database
    //!
    void NonInteractingOflModelData::ReadFromSql(
        const std::string tableName,    //!<    Name of SQL table where parameters are contained
        const std::string fileName,     //!<    Name of file where parameters are contained
        const iSize_t sqlId,            //!<    Sql identifier of the parameter set
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        utilities::Sqlite sql(fileName, _READ_EXISTING_);
        mpi.m_exitFlag = !sql.IsOpen();
        if(mpi.m_exitFlag)    
        {
            return;
        }
        utilities::SqliteRow sqlRow =  sql.RetrieveIdFromTable(tableName, sqlId);
        if(sqlRow.GetLength() == 0)
        {
            mpi.m_exitFlag = true;
            return;
        }
        else
	    {
	        utilities::cout.MainOutput()<<"\n\tReading model parameters from SQLite database "<<fileName<<std::endl;
	    }
        mpi.m_exitFlag = sqlRow.GetValue("Theta", m_theta); 
        if(mpi.m_exitFlag)    
        {
            return;
        }
        mpi.m_exitFlag = sqlRow.GetValue("V0", m_V0);
        if(mpi.m_exitFlag)    
        {
            return;
        }
        mpi.m_exitFlag = sqlRow.GetValue("Epsilon", m_epsilon);
        if(mpi.m_exitFlag)    
        {
            return;
        }
        mpi.m_exitFlag = sqlRow.GetValue("Kappa", m_kappa);
        if(mpi.m_exitFlag)    
        {
            return;
        }
        mpi.m_exitFlag = sqlRow.GetValue("Mass", m_mass);
        if(mpi.m_exitFlag)    
        {
            return;
        }
        return;
    }

    //!
    //! Default MPI function to synchronise all variables with those on the 
    //! given node
    //!
    void NonInteractingOflModelData::MpiSynchronize(
        const iSize_t nodeId,               //!<    Node id to sync with
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        mpi.Sync(&m_theta, 1, nodeId);
        mpi.Sync(&m_V0, 1, nodeId);
        mpi.Sync(&m_epsilon, 1, nodeId);
        mpi.Sync(&m_kappa, 1, nodeId);
        mpi.Sync(&m_mass, 1, nodeId);
    }
}   //  End namespace diagonalization
