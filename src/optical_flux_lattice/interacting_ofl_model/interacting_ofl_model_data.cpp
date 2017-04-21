////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                                                 
//!                                                                             
//!	 \file
//!     This file defines a class to store parameters defining an optical flux 
//!     lattice interacting model. See e.g. PRL 109, 265301 (2012)
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
#include "interacting_ofl_model_data.hpp"

namespace diagonalization
{
    //!
    //! Default Constructor
    //!
    InteractingOflModelData::InteractingOflModelData()
    :   m_dimX(2),
        m_dimY(2),
        m_nbrParticles(4),
        m_kxTot(0),
	    m_kyTot(0),
	    m_offsetX(0.0),
	    m_offsetY(0.0),
        m_outPath("./"),
        m_inPath("./"),
        m_outFileName("out"),
        m_nbrEigenvaluesToFind(4),
        m_interactionStrength(1.0),
        m_blockDiagonalize(false),
        m_useSql(false),
        m_sqlId(0),
        m_sqlName(""),
        m_sqlTableName(""),
        m_initialVectorFile("initVector.bin"),
        m_finalVectorFile("finalVector.bin"),
        m_tableFormat(_ARRAY_),
        m_setTableFormat(_ARRAY_),
        m_method(_FULL_),
        m_useWannierBasis(false),
	    m_termTablesBuilt(false),
	    m_magnetizationCalculated(false),
        m_fockBasisBuilt(false),
        m_hamiltonianBuilt(false),
        m_hamiltonianDiagonalized(false)
    {}

    //!
    //! Destructor
    //!
    InteractingOflModelData::~InteractingOflModelData()
    {}

    //!
    //! Copy constructor
    //!
    InteractingOflModelData::InteractingOflModelData(
        const InteractingOflModelData& other) //!<    Instance of the object
    :   m_dimX(other.m_dimX),
        m_dimY(other.m_dimY),
        m_nbrParticles(other.m_nbrParticles),
        m_kxTot(other.m_kxTot),
	    m_kyTot(other.m_kyTot),
	    m_offsetX(other.m_offsetX),
	    m_offsetY(other.m_offsetY),
        m_outPath(other.m_outPath),
        m_inPath(other.m_inPath),
        m_outFileName(other.m_outFileName),
        m_nbrEigenvaluesToFind(other.m_nbrEigenvaluesToFind),
        m_interactionStrength(other.m_interactionStrength),
        m_blockDiagonalize(other.m_blockDiagonalize),
        m_useSql(other.m_useSql),
        m_sqlId(other.m_sqlId),
        m_sqlName(other.m_sqlName),
        m_sqlTableName(other.m_sqlTableName),
        m_initialVectorFile(other.m_initialVectorFile),
        m_finalVectorFile(other.m_finalVectorFile),
        m_tableFormat(other.m_tableFormat),
        m_setTableFormat(other.m_setTableFormat),
        m_method(other.m_method),
        m_useWannierBasis(other.m_useWannierBasis),
	    m_termTablesBuilt(other.m_termTablesBuilt),
	    m_magnetizationCalculated(other.m_magnetizationCalculated),
        m_fockBasisBuilt(other.m_fockBasisBuilt),
        m_hamiltonianBuilt(other.m_hamiltonianBuilt),
        m_hamiltonianDiagonalized(other.m_hamiltonianDiagonalized)
    {}

    //!
    //! Overload assignment operator
    //!
    InteractingOflModelData& InteractingOflModelData::operator=(
        const InteractingOflModelData& other)   //!<    Instance of the object
    {
        m_dimX = other.m_dimX;
        m_dimY = other.m_dimY;
        m_offsetX = other.m_offsetX;
        m_offsetY = other.m_offsetY;
        m_nbrParticles = other.m_nbrParticles;
        m_kxTot = other.m_kxTot;
	    m_kyTot = other.m_kyTot;
        m_outPath = other.m_outPath;
        m_inPath = other.m_inPath;
        m_outFileName = other.m_outFileName;
        m_nbrEigenvaluesToFind = other.m_nbrEigenvaluesToFind;
        m_interactionStrength = other.m_interactionStrength;
        m_blockDiagonalize = other.m_blockDiagonalize;
        m_useSql = other.m_useSql;
        m_sqlId = other.m_sqlId;
        m_sqlName = other.m_sqlName;
        m_sqlTableName = other.m_sqlTableName;
        m_initialVectorFile = other.m_initialVectorFile;
        m_finalVectorFile = other.m_finalVectorFile;
        m_tableFormat = other.m_tableFormat;
        m_setTableFormat = other.m_setTableFormat;
        m_method = other.m_method;
        m_useWannierBasis = other.m_useWannierBasis;
        m_termTablesBuilt = other.m_termTablesBuilt;
        m_fockBasisBuilt = other.m_fockBasisBuilt;
        m_hamiltonianBuilt = other.m_hamiltonianBuilt;
        m_hamiltonianDiagonalized = other.m_hamiltonianDiagonalized;
        m_magnetizationCalculated = other.m_magnetizationCalculated;
        return *this;
    }

    //!
    //! MPI function to synchronise all variables with those on the given node
    //! 
    void InteractingOflModelData::MpiSynchronize(
        const int nodeId,                   //!<    Node id to sync with
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        mpi.Sync(m_inPath, nodeId);
        mpi.Sync(m_outPath, nodeId);
        mpi.Sync(m_outFileName, nodeId);
        mpi.Sync(m_sqlName, nodeId);
        mpi.Sync(m_sqlTableName, nodeId);
        mpi.Sync(m_initialVectorFile, nodeId);
        mpi.Sync(m_finalVectorFile, nodeId);
        mpi.Sync(&m_dimX, 1, nodeId);
        mpi.Sync(&m_dimY, 1, nodeId);
        mpi.Sync(&m_offsetX, 1, nodeId);
        mpi.Sync(&m_offsetY, 1, nodeId);
        mpi.Sync(&m_nbrParticles, 1, nodeId);
        mpi.Sync(&m_kxTot, 1, nodeId);
        mpi.Sync(&m_kyTot, 1, nodeId);
        mpi.Sync(&m_nbrEigenvaluesToFind, 1, nodeId);
        mpi.Sync(&m_interactionStrength, 1, nodeId);
        mpi.Sync(&m_useSql, 1, nodeId);
        mpi.Sync(&m_sqlId, 1, nodeId); 
        int useHash = m_setTableFormat;
        mpi.Sync(&useHash, 1, nodeId);
        if(useHash)
        {
            m_setTableFormat = _HASH_;
        }
        else
        {
            m_setTableFormat = _ARRAY_;
        }
        int method = 0 ;
        if(_FULL_ == m_method)
        {
            method = 0;
        }
        else if(_LANCZOS_ == m_method)
        {
            method = 1;
        }
        mpi.Sync(&method, 1, nodeId);
        if(0 == method)
        {
            m_method = _FULL_;
        }
        else if(1 == method)
        {
            m_method = _LANCZOS_;
        }
        mpi.Sync(&m_blockDiagonalize, 1, nodeId);
        mpi.Sync(&m_useWannierBasis, 1, nodeId);
        mpi.Sync(&m_termTablesBuilt, 1, nodeId);
        mpi.Sync(&m_fockBasisBuilt, 1, nodeId);
        mpi.Sync(&m_hamiltonianBuilt, 1, nodeId);
        mpi.Sync(&m_hamiltonianDiagonalized, 1, nodeId);
        mpi.Sync(&m_magnetizationCalculated, 1, nodeId); 
        if(m_useWannierBasis)
        {           
            m_setTableFormat = _HASH_;
        }
	    MPI_Barrier(mpi.m_comm);
    }
    
    //!
    //! Constructor for the InteractingOflModelData struct
    //! from command line arguments
    //!
    InteractingOflModelData::InteractingOflModelData(
        boost::program_options::variables_map* optionList,  
                                                    //!<    Command line arguments list
        utilities::MpiWrapper& mpi)                 //!<    Instance of the mpi wrapper class
        :
        m_blockDiagonalize(false),
        m_tableFormat(_ARRAY_), 
	    m_termTablesBuilt(false),
        m_hamiltonianBuilt(false),
        m_hamiltonianDiagonalized(false)
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        {
	            bool diagonalize;
	            GetOption(optionList, diagonalize, "diagonalize", _LINE_, mpi);
	            GetOption(optionList, m_dimX, "kx", _LINE_, mpi);
	            GetOption(optionList, m_dimY, "ky", _LINE_, mpi);
	            GetOption(optionList, m_offsetX, "kx-shift", _LINE_, mpi);
	            GetOption(optionList, m_offsetY, "ky-shift", _LINE_, mpi);
	            //  Note the offset values set from the command line will be overrided
                //  by the values optionally provided from an SQL database
	            GetOption(optionList, m_nbrParticles, "nbr", _LINE_, mpi);
	            GetOption(optionList, m_outPath, "out-path", _LINE_, mpi);
	            GetOption(optionList, m_inPath, "in-path", _LINE_, mpi);
	            GetOption(optionList, m_nbrEigenvaluesToFind, "nbr-eigenvalues", _LINE_, mpi);
	            GetOption(optionList, m_useSql, "use-sql", _LINE_, mpi);
	            GetOption(optionList, m_sqlName, "sql-name", _LINE_, mpi);
	            GetOption(optionList, m_sqlId, "sql-id", _LINE_, mpi);
	            GetOption(optionList, m_sqlTableName, "sql-table-name", _LINE_, mpi);
	            GetOption(optionList, m_initialVectorFile, "initial-file", _LINE_, mpi);
	            GetOption(optionList, m_finalVectorFile, "final-file", _LINE_, mpi);
	            bool useHash;
	            GetOption(optionList, useHash, "use-hash", _LINE_, mpi);
	            m_setTableFormat = myOptions::GetTermStorageType(useHash, mpi);
	            bool useFile;
                GetOption(optionList, useFile, "use-params-file", _LINE_, mpi);
	            if(mpi.m_exitFlag) 
	            {
	                goto escape;
	            }
                //  Need to read the mass parameter and interaction strength parameter
                //  from a text or sql file.
                std::stringstream fileName;
                fileName.str("");
                if(m_useSql)
                {
                    //  If sql option is set, then look for model data in an sqllite file
                    fileName << m_inPath << m_sqlName;
                    this->ReadFromSql(m_sqlTableName,fileName.str(), m_sqlId,diagonalize, mpi);
                    if(mpi.m_exitFlag) 
                    {
                        goto escape;
                    }
                }
                else if(useFile)
                {
                    std::string paramsFileName;
                    GetOption(optionList, paramsFileName, "params-file", _LINE_, mpi);
                    fileName << m_inPath << paramsFileName;
                    this->ReadFromFile(fileName.str(), mpi);
                    if(mpi.m_exitFlag) goto escape;
                    //  Generate a generic output file name  
                    fileName.str("");
                    fileName << "optical_flux_model_n_" << m_nbrParticles << "_kx_" << m_dimX << "_ky_" << m_dimY;
                    m_outFileName = fileName.str();
                }
                else
                {
                    double mass;
                    GetOption(optionList, mass, "mass", _LINE_, mpi);
                    GetOption(optionList, m_interactionStrength, "interaction-strength", _LINE_, mpi);
                    this->RescaleInteractionStrength(mass);   
                }
                //  Update the interaction strength with a factor 2*Pi/m as required
                utilities::cout.MainOutput()<<"\n\tMANY PARTICLE HAMILTONIAN PARAMETERS:"<<std::endl;
                utilities::cout.MainOutput()<<"\n\t\tdimX: \t\t\t"<<m_dimX<<std::endl;
                utilities::cout.MainOutput()<<"\t\tdimY: \t\t\t"<<m_dimY<<std::endl;
                utilities::cout.MainOutput()<<"\t\toffsetX: \t\t"<<m_offsetX<<std::endl;
                utilities::cout.MainOutput()<<"\t\toffsetY: \t\t"<<m_offsetY<<std::endl;
                utilities::cout.MainOutput()<<"\t\tnbrParticles: \t\t"<<m_nbrParticles<<std::endl;
                utilities::cout.MainOutput()<<"\t\tnbr eigs to find: \t"<<m_nbrEigenvaluesToFind<<std::endl;
                utilities::cout.MainOutput()<<"\t\tg_2D: \t\t\t"<<m_interactionStrength<<std::endl;
                utilities::cout.MainOutput()<<"\t\toutPath: \t\t"<<m_outPath<<std::endl;
                utilities::cout.MainOutput()<<"\t\tinPath: \t\t"<<m_inPath<<std::endl;
                if(m_useSql)
                {
                    utilities::cout.SecondaryOutput()<<"\t\tUsing Sqlite file \t"<<m_sqlName<<std::endl;
                    utilities::cout.SecondaryOutput()<<"\t\tSqlite parameters ID\t"<<m_sqlId<<std::endl;
                }
                else
                {
                    utilities::cout.SecondaryOutput()<<"\t\tOut file name \t\t"<<m_outFileName<<std::endl;
                }
                int methodCode;
                GetOption(optionList, methodCode, "method", _LINE_, mpi);
                m_method = myOptions::GetDiagonalizationMethod(methodCode, mpi);
                if(mpi.m_exitFlag) 
                {
                    goto escape;
                }
                if(_FULL_ == m_method)
                {
                    utilities::cout.SecondaryOutput()<<"\n\tDiagonalization method: Full (LAPACK)"<<std::endl;
                }
                else if(_LANCZOS_ == m_method)
                {
                    utilities::cout.SecondaryOutput()<<"\n\tDiagonalization method: Lanczos (ARPACK)"<<std::endl;
                }
                int basisCode;
                GetOption(optionList, basisCode, "basis" , _LINE_, mpi);
                m_useWannierBasis = myOptions::GetBasisType(basisCode, mpi);
                if(mpi.m_exitFlag) 
                {
                    goto escape;
                }
                if(m_useWannierBasis)
                {           
                    utilities::cout.SecondaryOutput()<<"\n\tUsing Hybrid Localized Wannier basis (forces hash type term storage)"<<std::endl;
                    m_setTableFormat = _HASH_;
                }
            }
            escape:
            {}
        }
	    mpi.ExitFlagTest();
	    this->MpiSynchronize(0, mpi);
    }

    //!
    //! Function to read parameters from a file
    //!
    void InteractingOflModelData::ReadFromFile(
        const std::string fileName,     //!<    Name of file where parameters are contained
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        if(mpi.m_id == 0)
        {
            std::ifstream f_param;
            io::fileFormat_t format = io::_TEXT_;
            utilities::GenFileStream(f_param, fileName, format, mpi);
            if(mpi.m_exitFlag)
            {
                return;
            }
            std::string line;
            while(getline(f_param, line))
            {
                if(line[0]!='#')	break;
            }
            for(int i=0; i<4; ++i)
            {   
                //  Discard first 4 lines
                getline(f_param, line);
            }
            double mass;
            f_param >> mass >> m_interactionStrength;
            f_param.close();
            this->RescaleInteractionStrength(mass);
        }
        return;
    }

    //!
    //! Function to read parameters from a sql database
    //!
    void InteractingOflModelData::ReadFromSql(
        const std::string tableName,    //!<    Name of SQL table where parameters are contained
        const std::string fileName,     //!<    Name of file where parameters are contained
        const iSize_t sqlId,            //!<    Sql identifier of the parameter set
        const bool diagonalize,         //!<    Option to diagonalize the Hamiltonian or not
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        utilities::Sqlite sql(fileName, sql::_READ_EXISTING_);
        mpi.m_exitFlag = !sql.IsOpen();
        if(mpi.m_exitFlag)    
        {
            return;
        }
        utilities::SqliteRow sqlRow = sql.RetrieveIdFromTable(tableName, sqlId); 
        if(sqlRow.GetLength() == 0)
        {
            mpi.m_exitFlag = true;
            return;
        }
        //  Check if the calculated has already been completed (as recorded by
        //  the "completed" column in the SQL table
        int isAlreadyCalculated = 0;
        mpi.m_exitFlag = sqlRow.GetValue("Completed", isAlreadyCalculated);
        if(mpi.m_exitFlag)    
        {
            return;
        }
        if(isAlreadyCalculated && diagonalize)
        {
            std::cerr<<"Calculated already completed (according to SQL table). SKIPPING"<<std::endl;
            mpi.m_exitFlag = true;
            return;
        }
        if(!isAlreadyCalculated && !diagonalize)
        {
            std::cerr<<"ERROR Requested Spectrum not calculated (according to SQL table)."<<std::endl;
            mpi.m_exitFlag = true;
            return;
        }
	    double mass = 1.0;  //  set default value
	    mpi.m_exitFlag = sqlRow.GetValue("Mass", mass);	
        if(mpi.m_exitFlag)    
        {
            return;
        }   
	    mpi.m_exitFlag = sqlRow.GetValue("Interaction", m_interactionStrength);
        if(mpi.m_exitFlag)    
        {
            return;
        }
	    mpi.m_exitFlag = sqlRow.GetValue("OutFileName", m_outFileName);
        if(mpi.m_exitFlag)    
        {
            return;
        }
        mpi.m_exitFlag = sqlRow.GetValue("OffsetX", m_offsetX);	
        if(mpi.m_exitFlag)    
        {
            return;
        }
        mpi.m_exitFlag = sqlRow.GetValue("OffsetY", m_offsetY);	
        if(mpi.m_exitFlag)    
        {
            return;
        }
        this->RescaleInteractionStrength(mass);
        return;
    }

    //!
    //! Rescale the interaction strength by mass/2pi
    //!
    void InteractingOflModelData::RescaleInteractionStrength(
        const double mass)  //!<    Particle mass in the model
    {
        m_interactionStrength *= 2.0*PI/mass;
    }

    //!
    //! Generate an output file name base, using struct parameters
    //!
    std::string InteractingOflModelData::MakeBaseFileName(
        const io::io_t io)  //!<    Specify the creation of an input or output file name
        const
    {
        std::stringstream fileName;   
        fileName.str("");
        if(io::_OUT_ == io)
        {
            fileName<<m_outPath;
        }
        else if(io::_IN_ == io)
        {
            fileName<<m_inPath;
        }
        fileName<<"/"<<m_outFileName;
        if(m_blockDiagonalize)
        {
            if(m_useWannierBasis)
            {
                fileName<<"_sector_"<<m_kyTot;
            }
            else
            {
                fileName<<"_sector_"<<m_kxTot<<"_"<<m_kyTot;
            }
        }
        return fileName.str();
    }

    //!
    //! Generate an output file name base, using struct parameters
    //!
    void InteractingOflModelData::UpdateSqlStatus(
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        //  Update SQL table to indicate that the result 
        //  has been calculated and update the time entered field
        if(0 == mpi.m_id && m_useSql)    //  FOR THE MASTER NODE
        {
            std::stringstream fileName;
            fileName.str("");
            fileName<<m_inPath<<m_sqlName;
            utilities::Sqlite sql(fileName.str(), sql::_EDIT_EXISTING_);
            {
                utilities::SqliteVariable var;
                var.SetValues("Completed", (int)1);
                sql.UpdateTableEntry(m_sqlTableName, m_sqlId, var);
            }
            //  Record the PBS job id, if known
            if(getenv("PBS_JOBID")!=NULL)
            {
                utilities::SqliteVariable var;
                var.SetValues("jobId",(int)atoi(getenv("PBS_JOBID")));
                sql.UpdateTableEntry(m_sqlTableName, m_sqlId, var);
            }
            //  Record the SLURM job id, if known
		    if(getenv("SLURM_JOB_ID")!=NULL)
            {
                utilities::SqliteVariable var;
                var.SetValues("jobId", (int)atoi(getenv("SLURM_JOB_ID")));
                sql.UpdateTableEntry(m_sqlTableName, m_sqlId, var);
            }
            sql.CopyTableToTextFile(m_sqlTableName, fileName.str());
        }
    }
}   //  End namespace diagonalization
