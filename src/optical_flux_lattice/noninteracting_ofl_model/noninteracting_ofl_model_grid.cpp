////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store solutions to a non-interacting 
//!     optical flux lattice model at different k-space points.
//!     See e.g. PRL 109, 265301 (2012).
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
#include "noninteracting_ofl_model_grid.hpp"

namespace diagonalization
{
    //!
    //! This function combines a set of 2D indexes into a single index
    //!
    kState_t NonInteractingOflModelGrid::CombineIndexes(
        const kState_t x,           //!<    x index
        const kState_t y,           //!<    y index
        const iSize_t dim)          //!<    Number of points on the grid
        const
    {
        return x*dim + y;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief  Generate a table single-particle Hamiltonian solutions for a given 
    //! k-space grid covering the Brillouin zone with dimX by dimY points.
    //!
    //! Note: pass the address of an array of NonInteractingOflModel classes
    //! as the Bloch wave function look-up table
    //!
    //! Note: This function is parallelized
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModelGrid::GenerateBlochWaveFunctionTable(
        const iSize_t cutOff,                   //!<    Cut-off in range of reciprocal lattice vectors
        const iSize_t dimX,                     //!<    Number of k-space points in G1-direction
        const iSize_t dimY,                     //!<    Number of k-space points in G2-direction
        const double offsetX,                   //!<    X offset in the map of k-space
        const double offsetY,                   //!<    Y offset in the map of k-space
        NonInteractingOflModel* blochTable,     //!<    Look-up table of Bloch wave functions  
        const NonInteractingOflModelData* params, //!<    Parameter set defining the model
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING BLOCH WAVE FUNCTION TABLE (cut-off range "<<cutOff<<")"<<std::endl<<std::endl;
        }
        //  Generate a list of NonInteractingOflModel objects to store
        //  tables of Bloch wave function data
        mpi.DivideTasks(mpi.m_id,dimX*dimY, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
        NonInteractingOflModel oflModel(*params);
        oflModel.m_xBandCutOff = cutOff;
        oflModel.m_yBandCutOff = cutOff;
        oflModel.m_nbrBands = 1;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.DebuggingInfo()<<"\n\t- ALLOCATING ~"<<dimX*dimY*oflModel.GetDifferenceMatrixDimension()*sizeof(dcmplx)/(1024.0*1024.0)<<" mb";
            utilities::cout.DebuggingInfo()<<" to store Bloch wave function loop-up table.\n"<<std::endl;
        }
        oflModel.PopulateDifferenceMatrix();
        utilities::LatticeVector2D<double> kVector(offsetX, offsetY);
        utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1, (double)1/dimX);
        utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2, (double)1/dimY);
        iSize_t counter  = 0;
        utilities::LoadBar progress;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            progress.Initialize(mpi.m_lastTask - mpi.m_firstTask+1);
        }        
        //  Allocate memory to store part of the result on each node
        std::vector<NonInteractingOflModel> temp;
        NonInteractingOflModel dummy;
        temp.reserve(ceil((double)mpi.m_nbrProcs/(dimX*dimY)));
        MPI_Barrier(mpi.m_comm);
        // Parallelize the the for loop by dividing the k-space grid calculations
        // up over a number of processors. 
        for(iSize_t x=0; x<dimX; ++x, kVector-=reciprocalLattice::G2, kVector += increment1)   //  Map out G1 direction
        {
            for(iSize_t y=0; y<dimY; ++y, kVector += increment2, ++counter)   //  Map out G2 direction
            {
                if(mpi.m_firstTask <= (int)counter && mpi.m_lastTask >= (int)counter)
                {
                    oflModel.UpdateDifferenceDiagonal(kVector);
                    oflModel.DiagonalizeDifferenceMatrix();
                    //  Set the m_differenceMatrixPopulated to false so that the
                    //  full difference matrix is not coppied
                    oflModel.m_differenceMatrixPopulated = false;
                    temp.push_back(oflModel);
                    //  Reset the flag so that the difference matrix can be updated
                    oflModel.m_differenceMatrixPopulated = true;
                }
                else
                {
                    //  Append empty class to the list to assist with the MPI sync function below
                    temp.push_back(dummy);
                }
                if(0 == mpi.m_id && mpi.m_lastTask >= (int)counter)	// FOR THE MASTER NODE
	            {
                    progress.Display(counter+1);
                }
            }
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<std::endl;
        }
        //  Gather up all the calculated values and put them in a copy of the Bloch
        //  table on each node. Example (for objects T0,T1,T2,T3,T4) 
        //
        //          NODE0       |       NODE1    |       NODE2
        //
        //  INPUT:  T0 T1 _ _ _ |  _ _ T2 T3 _   |  _ _ _ _ T4 
        //  (temp array)
        //
        //  OUTPUT: all T0 to T4|  all T0 to T4  |  all T0 to T4
        //  (blochTable)
        //
        MPI_Barrier(mpi.m_comm);
        counter = 0;
        int proc_remain = std::min(mpi.m_nbrProcs, (int)(dimX*dimY));
        int taskRemain = dimX*dimY;
        NonInteractingOflModel* p_table = blochTable;
        int taskCounter = 0;
        for(int node=0; node<mpi.m_nbrProcs; ++node)
        {
            int taskProc = (int)floor((double)taskRemain/proc_remain+0.5);
            proc_remain--;
            taskRemain -= taskProc;
            for(int p=0; p<taskProc; ++p, ++taskCounter, ++p_table)
            {
                p_table->MpiSynchrnoizedCopy(temp[taskCounter], node, mpi);
            } 
        }
        return;
    }
    
    //!
    //! Numerically calculates the non-interacting band structure of the 
    //! Hamiltonian given in e.g. eq (1) of PRL 109,265301 in a single Brillouin zone
    //!
    void NonInteractingOflModelGrid::CalculateBandstructure(
        double* lowerBand,          //!<    List of energy eigenvalues in the lowest band
                                    //!     (returned on node 0 only)
        double* secondBand,         //!<    List of energy eigenvalues in the second band
                                    //!     (returned on node 0 only)
        double& bandWidth,          //!<    Return the difference between the maximum
                                    //!     and minimum energy in the lowest band
                                    //!     (returned on node 0 only)
        double& bandGap,            //!<    Return the minimum difference in energy
                                    //!     between the lowest and second lowest band
                                    //!     (returned on node 0 only)
        const bool getEnergyData,   //!<    Set to true to populate the arrays in 
                                    //!     lowerBand,secondBand,bandWidth,bandGap
                                    //!     (returned on node 0 only)
        double* magnetization,      //!<    Return a list of magnetization values
                                    //!     (returned on node 0 only)
        const bool getMagnetizationData,
                                    //!<    Set to true to populate the magnatization array
        const iSize_t dimX,         //!<    Number of k-space points in G1-direction
        const iSize_t dimY,         //!<    Number of k-space points in G2-direction
        const double offsetX,       //!<    X offset in the map of k-space
        const double offsetY,       //!<    Y offset in the map of k-space
        NonInteractingOflModel* currentModel,
                                    //!<    Class containing the model data
        utilities::MpiWrapper& mpi) //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.SecondaryOutput()<<"\n\t- GENERATING BAND STRUCTURE FOR CURRENT MODEL PARAMETERS"<<std::endl;
        }
        currentModel->PopulateDifferenceMatrix();
        utilities::LatticeVector2D<double> kVector(offsetX,offsetY);
        utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1,(double)1/dimX);
        utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2,(double)1/dimY);
        mpi.DivideTasks(mpi.m_id,dimX*dimY,mpi.m_nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);
        iSize_t counter  = 0;
        utilities::LoadBar progress;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            progress.Initialize(mpi.m_lastTask - mpi.m_firstTask+1);
        }
        for(iSize_t x=0; x<dimX; ++x, kVector-=reciprocalLattice::G2, kVector += increment1)   //  Map out G1 direction
        {
            for(iSize_t y=0; y<dimY; ++y, kVector += increment2, ++counter)   //  Map out G2 direction
            {
                if(mpi.m_firstTask <= (int)counter && mpi.m_lastTask >= (int)counter)
                {      
                    currentModel->UpdateDifferenceDiagonal(kVector);
                    currentModel->DiagonalizeDifferenceMatrix();
                    if(getEnergyData)
                    {
                        lowerBand[x*dimY+y]  = currentModel->m_eigenvalues[0];
                        secondBand[x*dimY+y] = currentModel->m_eigenvalues[1];
                    }
                    if(getMagnetizationData)
                    {
                        magnetization[x*dimY+y] = currentModel->ApplySigmaZ(0);
                    }
                }
                if(0 == mpi.m_id && mpi.m_lastTask >= (int)counter)	// FOR THE MASTER NODE
	            {
                    progress.Display(counter+1);
                }
            }
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<std::endl;
        }
        //  Call MPI Gather on the energy and magnetization arrays
        double* gatherBuffer = 0;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        gatherBuffer = new double[dimX*dimY];
	    }
        int bufferSize = mpi.m_lastTask - mpi.m_firstTask +1;
        if(getEnergyData)
        {
            MPI_Status status;
            mpi.Gather(lowerBand+mpi.m_firstTask, bufferSize, gatherBuffer, dimX*dimY, 0, mpi.m_comm, status);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                memcpy(lowerBand, gatherBuffer, dimX*dimY*sizeof(double));
            }
            mpi.Gather(secondBand+mpi.m_firstTask, bufferSize, gatherBuffer, dimX*dimY, 0, mpi.m_comm, status);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                memcpy(secondBand, gatherBuffer, dimX*dimY*sizeof(double));
            }
        }
        if(getMagnetizationData)
        {
            MPI_Status status;
            mpi.Gather(magnetization+mpi.m_firstTask, bufferSize, gatherBuffer, dimX*dimY, 0, mpi.m_comm, status);
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                memcpy(magnetization, gatherBuffer, dimX*dimY*sizeof(double));
            }
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        delete[] gatherBuffer;
	    }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            if(getEnergyData)
            {
                //  Calculate band width of lowest band    
                bandWidth = *std::max_element(lowerBand,lowerBand+dimX*dimY)
                          - *std::min_element(lowerBand,lowerBand+dimX*dimY);
                double bandGapList[dimX*dimY];
                for(iSize_t i=0;i<dimX*dimY;++i)
                {
                    bandGapList[i] = secondBand[i] - lowerBand[i];   
                }
                bandGap = *std::min_element(bandGapList, bandGapList+dimX*dimY);
            }
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Plots the single particle band structure of the 
    //! Hamiltonian given in e.g. eq (1) of PRL 109,265301 for a given sample
    //! of k-space.
    //!
    //! Note: This function is parallelized
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModelGrid::PlotBandstructure(
        boost::program_options::variables_map* optionList,     
                                                    //!<    Command line arguments
        utilities::MpiWrapper& mpi)                 //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.MainOutput()<<"\n\tGENERATING SINGLE PARTICLE BAND STRUCTURE DATA..."<<std::endl;
        }
        iSize_t dimX;
        iSize_t dimY;
        double offsetX;
        double offsetY;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        GetOption(optionList, dimX, "x-grid", _LINE_, mpi);
	        GetOption(optionList, dimY, "y-grid", _LINE_, mpi);
	        GetOption(optionList, offsetX, "x-grid-shift", _LINE_, mpi);
	        GetOption(optionList, offsetY, "x-grid-shift", _LINE_, mpi);
        }
        mpi.Sync(&dimX, 1, 0);
        mpi.Sync(&dimY, 1, 0);
        mpi.Sync(&offsetX, 1, 0);
        mpi.Sync(&offsetY, 1, 0);
        NonInteractingOflModel currentModel(optionList, mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.SecondaryOutput()<<"\n\tGRID: "<<dimX<<" x "<<dimY<<std::endl;
            utilities::cout.SecondaryOutput()<<"\n\tG_X CUT-OFF = "<<currentModel.m_xBandCutOff<<std::endl;
            utilities::cout.SecondaryOutput()<<"\n\tG_Y CUT-OFF = "<<currentModel.m_yBandCutOff<<std::endl;
        }
        if(2 != currentModel.m_nbrBands)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                std::cerr<<"WARNING: ONLY GENERATING DATA FOR THE LOWEST TWO BANDS"<<std::endl;
            }
            currentModel.m_nbrBands = 2;
        }
        //  Allocate memory to store band energy values
        double lowerBand[dimX*dimY];
        double secondBand[dimX*dimY];
        double bandGap   = 0.0;
        double bandWidth = 0.0;
        double dummy;
        this->CalculateBandstructure(lowerBand, secondBand, bandWidth, bandGap, true, &dummy, 
                                     false, dimX, dimY, offsetX, offsetY, &currentModel, mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            std::stringstream dataFileName;
            dataFileName.str("");
            std::string outPath;
            GetOption(optionList, outPath, "out-path", _LINE_, mpi);
            dataFileName << optionList << "/band_structure_kx_" << dimX << "_ky_" << dimY;
            bool useSql;
            GetOption(optionList, useSql, "use-sql", _LINE_, mpi);
            iSize_t sqlId;
            GetOption(optionList, sqlId, "sql-id", _LINE_, mpi);
            if(useSql)
            {
                dataFileName << "_id_" << sqlId;
            }
            dataFileName << ".dat";
            std::ofstream f_tmp;
            f_tmp.open(dataFileName.str().c_str(), std::ios::out);
            double* p_lower = lowerBand;
            double* p_second = secondBand;
            utilities::LatticeVector2D<double> kVector(offsetX, offsetY);
            utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1, (double)1/dimX);
            utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2, (double)1/dimY); 
            f_tmp<<dimX<<"\n";
            f_tmp<<dimY<<"\n";
            f_tmp<<std::setw(15)<<std::left<<"kx values"<<"\t"<<std::setw(15)<<std::left<<"ky values"<<"\t"<<std::setw(20)<<std::left<<"0th band energy"<<"\t"<<std::setw(20)<<std::left<<"1st band energy"<<"\n";
            for(iSize_t x=0; x<dimX; ++x, kVector-=reciprocalLattice::G2, kVector += increment1)
            {
                for(iSize_t y=0; y<dimY; ++y, ++p_lower, ++p_second, kVector += increment2)
                {
                    f_tmp<<std::setw(15)<<std::setprecision(10)<<std::left<<kVector.GetKx()<<"\t"<<std::setw(15)<<std::left<<kVector.GetKy()<<"\t"<<std::setw(20)<<std::left<<std::setprecision(15)<<*(p_lower)<<"\t"<<std::setw(20)<<std::left<<std::setprecision(15)<<*(p_second)<<"\n";
                }
            }
            f_tmp.close();
            utilities::cout.SecondaryOutput()<<"\n\tData output to file "<<dataFileName.str()<<std::endl;
            bool plotBandstructure;
            GetOption(optionList, plotBandstructure, "plot-bandstructure", _LINE_, mpi);
            if(plotBandstructure)
            {
                utilities::cout.AdditionalInfo()<<"\n\tMAKING PLOT..."<<std::endl;
                std::stringstream plotFileName;
                plotFileName.str("");
                std::string outPath;
                GetOption(optionList, outPath, "out-path", _LINE_, mpi);
                plotFileName<<outPath<<"/band_structure_kx_"<<dimX<<"_ky_"<<dimY<<".pdf";
                std::stringstream pythonScript;
                pythonScript.str();
                pythonScript<<"#! //usr//bin//env python"<<_PYTHON_VERSION_<<"\n"
                "import matplotlib                              \n"
                "import numpy as np                             \n"
                "import math                                    \n"
                "matplotlib.use('Agg')                        \n\n"
                "import matplotlib.pyplot as plt              \n\n"
                "from mpl_toolkits.mplot3d import Axes3D      \n\n"
                "xData = []                                     \n"
                "yData = []                                     \n"   
                "lowerBand = []                                 \n"
                "secondBand = []                              \n\n"  
                "fin=open(\""<<dataFileName.str()<<"\")       \n\n"
                "fin.readline()                                 \n"
                "fin.readline()                                 \n"
                "fin.readline()                                 \n"
                "for line in fin:                               \n"
                "    line = line.split()                        \n"
                "    xData.append(float(line[0]))               \n"
                "    yData.append(float(line[1]))               \n"
                "    lowerBand.append(float(line[2]))           \n"
                "    secondBand.append(float(line[3]))          \n"
                "fin.close()                                  \n\n"
                "fig = plt.figure()                           \n\n"
                "ax = fig.gca(projection='3d')                \n\n"
                "ax.plot_trisurf(xData,yData,lowerBand)         \n"
                "ax.plot_trisurf(xData,yData,secondBand,color='r')\n"
                "plt.title(\"Band gap: "<<bandGap<<" ; Band width: "<<bandWidth<<"\")\n"
                "plt.savefig(\""<<plotFileName.str()<<"\") \n"
                "plt.close() \n";
                utilities::Script myScript;
                myScript.SetScript(pythonScript.str());
                myScript.Execute();
                utilities::cout.MainOutput()<<"\n\t...PLOT SUCESSFULLY GENERATED"<<std::endl;
                utilities::cout.MainOutput()<<"\n\tFind it at "<<plotFileName.str()<<std::endl;
            }
        }
        return;
    }

    //!
    //! Make a plot of the band width as a function of the lattice depth V0
    //!
    void NonInteractingOflModelGrid::PlotBandWidth(
        boost::program_options::variables_map* optionList,     
                                            //!<    Command line arguments
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.MainOutput()<<"\n\tGENERATING BAND WIDTH AND GAP VS V0/E_R DATA..."<<std::endl;
        }
        iSize_t dimX;
        iSize_t dimY;
        double offsetX;
        double offsetY;
        std::stringstream dataFileName;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        GetOption(optionList, dimX, "x-grid", _LINE_, mpi);
	        GetOption(optionList, dimY, "y-grid", _LINE_, mpi);
	        bool useSql;
	        GetOption(optionList, useSql, "use-sql", _LINE_, mpi);
	        std::string outPath;
	        GetOption(optionList, outPath, "out-path", _LINE_, mpi);
	        if(useSql)
            {
                std::stringstream fileName;
	            fileName.str("");
                //  If sql option is set, then look for model data in an sqlite file
                std::string inPath;
                std::string sqlName;
                std::string sqlTableName;
                GetOption(optionList, inPath, "in-path", _LINE_, mpi);
                GetOption(optionList, sqlName, "sql-name", _LINE_, mpi);
                GetOption(optionList, sqlTableName, "sql-table-name", _LINE_, mpi);
                fileName << inPath << sqlName;
                utilities::Sqlite sql(fileName.str(), sql::_READ_EXISTING_);
                if(!mpi.m_exitFlag)
                {
                    mpi.m_exitFlag = !sql.IsOpen();
                }
                iSize_t sqlId;
                GetOption(optionList, sqlId, "sql-id", _LINE_, mpi);
                utilities::SqliteRow sqlRow = sql.RetrieveIdFromTable(sqlTableName, sqlId); 
                if(sqlRow.GetLength() == 0)
                {
                    mpi.m_exitFlag = true;
                }
                //  Check if the calculated has already been completed (as recorded by
                //  the "completed" column in the SQL table
                int isAlreadyCalculated = 0;
                mpi.m_exitFlag = sqlRow.GetValue("Completed", isAlreadyCalculated);
                if(isAlreadyCalculated)
                {
                    std::cerr<<"Calculated already completed (according to SQL table). SKIPPING"<<std::endl;
                    mpi.m_exitFlag = true;
                }
                std::string outFileName;
                mpi.m_exitFlag = sqlRow.GetValue("OutFileName",outFileName);
                mpi.m_exitFlag = sqlRow.GetValue("OffsetX",offsetX);	
                mpi.m_exitFlag = sqlRow.GetValue("OffsetY",offsetY);	
                dataFileName.str("");
                dataFileName << outPath << outFileName << ".dat";
            }
            else
            {
                GetOption(optionList, offsetX, "x-grid-shift", _LINE_, mpi);
                GetOption(optionList, offsetY, "y-grid-shift", _LINE_, mpi);
                dataFileName.str("");
                dataFileName << outPath << "/band_width_kx_" << dimX << "_ky_" << dimY << ".dat";
            }
        }
        mpi.ExitFlagTest();
        mpi.Sync(&dimX, 1, 0);
        mpi.Sync(&dimY, 1, 0);
        mpi.Sync(&offsetX, 1, 0);
        mpi.Sync(&offsetY, 1, 0);
        NonInteractingOflModel currentModel(optionList, mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.SecondaryOutput()<<"\n\tGRID: "<<dimX<<" x "<<dimY<<std::endl;
            utilities::cout.SecondaryOutput()<<"\n\tG_X CUT-OFF = "<<currentModel.m_xBandCutOff<<std::endl;
            utilities::cout.SecondaryOutput()<<"\n\tG_Y CUT-OFF = "<<currentModel.m_yBandCutOff<<std::endl;
        }
        //  Allocate memory to store band energy values
        double lowerBand[dimX*dimY];
        double secondBand[dimX*dimY];
        const double V0Max = 5.0;
        const double V0grid = 0.1;
        const iSize_t gridMax = floor(V0Max/V0grid);
        //  Update the V0 parameter to the first V0 value
        currentModel.m_params.m_V0 = V0grid;
        //  Output the plot data as a .dat file
        std::ofstream f_tmp;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            f_tmp.open(dataFileName.str().c_str(), std::ios::out);
            f_tmp<<std::setw(15)<<std::left<<"V0/E_R values"<<"\t"<<std::setw(15)<<std::left<<"band width/E_R"<<"\t"<<std::setw(15)<<std::left<<"band gap/E_R"<<"\n";
        }
        //  Generate the band structure for each V0 value
        for(iSize_t v=0; v<=gridMax; ++v)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                utilities::cout.SecondaryOutput()<<"\n\tV0 = "<<currentModel.m_params.m_V0<<" (going up to "<<V0Max+V0grid<<")"<<std::endl;
            }
            double bandGap;
            double bandWidth;
            double dummy;
            this->CalculateBandstructure(lowerBand, secondBand, bandWidth, bandGap, true, &dummy, 
                                         false, dimX, dimY, offsetX, offsetY, &currentModel, mpi);
            
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
                f_tmp<<std::setw(15)<<std::left<<currentModel.m_params.m_V0/currentModel.GetRecoilEnergy()<<"\t"<<std::setw(15)<<std::left<<bandWidth/currentModel.GetRecoilEnergy()<<"\t"<<std::setw(15)<<std::left<<bandGap/currentModel.GetRecoilEnergy()<<"\n";
            }
            currentModel.m_params.m_V0 += V0grid;
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            f_tmp.close();
            utilities::cout.SecondaryOutput()<<"\n\tData output to file "<<dataFileName.str()<<std::endl;
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        bool plotBandWidth;
	        GetOption(optionList, plotBandWidth, "plot-band-width", _LINE_, mpi);
	        std::string outPath;
	        GetOption(optionList, outPath, "out-path", _LINE_, mpi);
            if(plotBandWidth)
            {
                utilities::cout.SecondaryOutput()<<"\n\tMAKING PLOT..."<<std::endl;
                std::stringstream plotFileName;
                plotFileName.str("");
                plotFileName<<outPath<<"/band_width_kx_"<<dimX<<"_ky_"<<dimY<<".pdf";
                std::stringstream pythonScript;
                pythonScript.str();
                pythonScript<<"#! //usr//bin//env python"<<_PYTHON_VERSION_<<"\n"
                "import matplotlib                              \n"
                "import numpy as np                             \n"
                "import math                                    \n"
                "matplotlib.use('Agg')                        \n\n"
                "import matplotlib.pyplot as plt              \n\n"
                "from mpl_toolkits.mplot3d import Axes3D      \n\n"
                "xData = []                                     \n"
                "bandWidth = []                                 \n"
                "bandGap = []                                 \n\n"  
                "fin=open(\""<<dataFileName.str()<<"\")       \n\n"
                "fin.readline()                                 \n"
                "for line in fin:                               \n"
                "    line = line.split()                        \n"
                "    xData.append(float(line[0]))               \n"
                "    bandWidth.append(float(line[1]))           \n"
                "    bandGap.append(float(line[2]))             \n"
                "fin.close()                                  \n\n"
                "p1, = plt.plot(xData,bandWidth)                \n"
                "p2, = plt.plot(xData,bandGap,c='r')            \n"
                "plt.title(\"Band width and gap vs $V_0/E_R$\") \n"
                "plt.legend([p1,p2],[\"Band Width/E_R\",\"Band Gap/E_R\"])\n"
                "plt.savefig(\""<<plotFileName.str()<<"\") \n"
                "plt.close() \n";
                utilities::Script myScript;
                myScript.SetScript(pythonScript.str());
                myScript.Execute();
                utilities::cout.MainOutput()<<"\n\t...PLOT SUCESSFULLY GENERATED"<<std::endl;   
                utilities::cout.MainOutput()<<"\n\tFind it at "<<plotFileName.str()<<std::endl;
            }
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        bool useSql;
	        GetOption(optionList, useSql, "use-sql", _LINE_, mpi);
	        std::string inPath;
	        GetOption(optionList, inPath, "in-path", _LINE_, mpi);
	        if(useSql)
            {
                std::stringstream fileName;
	            fileName.str("");
                //  If sql option is set, then look for model data in an sqlite file
                std::string sqlName;
                GetOption(optionList, sqlName, "sql-name", _LINE_, mpi);
                fileName << inPath << sqlName;
                std::string sqlTableName;
                iSize_t sqlId;
                GetOption(optionList, sqlTableName, "sql-table-name", _LINE_, mpi);
                GetOption(optionList, sqlId, "sql-id", _LINE_, mpi);
                utilities::Sqlite sql(fileName.str(), sql::_EDIT_EXISTING_);
                {
                    utilities::SqliteVariable var;
                    var.SetValues("Completed", (int)1);
                    sql.UpdateTableEntry(sqlTableName, sqlId, var);
                }
                //  Record the PBS job id, if known
                if(getenv("PBS_JOBID")!=NULL)
                {
                    utilities::SqliteVariable var;
                    var.SetValues("jobId",(int)atoi(getenv("PBS_JOBID")));
                    sql.UpdateTableEntry(sqlTableName, sqlId, var);
                }
                //  Record the SLURM job id, if known
		        if(getenv("SLURM_JOB_ID")!=NULL)
                {
                    utilities::SqliteVariable var;
                    var.SetValues("jobId", (int)atoi(getenv("SLURM_JOB_ID")));
                    sql.UpdateTableEntry(sqlTableName, sqlId, var);
                }
                //  Update the text copy of the table
                sql.CopyTableToTextFile(sqlTableName, fileName.str());
            }
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Make a plot of the magnetization for a given 2D k-space grid
    //!
    //! Note: This function is parallelized
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModelGrid::PlotMagnetization(
        boost::program_options::variables_map* optionList,  
                                                //!<    Command line arguments
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.MainOutput()<<"\n\tGENERATING SINGLE PARTICLE MAGNETIZATION DATA..."<<std::endl;
        }
        iSize_t dimX;
        iSize_t dimY;
        double offsetX;
        double offsetY;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        GetOption(optionList, dimX, "x-grid", _LINE_, mpi);
	        GetOption(optionList, dimY, "y-grid", _LINE_, mpi);
	        GetOption(optionList, offsetX, "x-grid-shift", _LINE_, mpi);
	        GetOption(optionList, offsetY, "y-grid-shift", _LINE_, mpi);
        }
        mpi.Sync(&dimX, 1, 0);
        mpi.Sync(&dimY, 1, 0);
        mpi.Sync(&offsetX, 1, 0);
        mpi.Sync(&offsetY, 1, 0);
        NonInteractingOflModel currentModel(optionList, mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.SecondaryOutput()<<"\n\tGRID: "<<dimX<<" x "<<dimY<<std::endl;
            utilities::cout.SecondaryOutput()<<"\n\tG_X CUT-OFF = "<<currentModel.m_xBandCutOff<<std::endl;
            utilities::cout.SecondaryOutput()<<"\n\tG_Y CUT-OFF = "<<currentModel.m_yBandCutOff<<std::endl;
        }
        //  Get magnetization data
        double dummy;
        double magnetization[dimX*dimY];
        this->CalculateBandstructure(&dummy, &dummy, dummy, dummy, false, magnetization, true, 
                                     dimX, dimY, offsetX, offsetY, &currentModel, mpi);
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            std::stringstream dataFileName;
            dataFileName.str("");
            std::string outPath;
            GetOption(optionList, outPath, "out-path", _LINE_, mpi);
            dataFileName << outPath << "/magnetization_map_kx_" << dimX << "_ky_" << dimY;
            bool useSql;
            GetOption(optionList, useSql, "use-sql", _LINE_, mpi);
            iSize_t sqlId;
            GetOption(optionList, sqlId, "sql-id", _LINE_, mpi);
            if(useSql)
            {
                dataFileName << "_id_" << sqlId;
            }
            dataFileName << ".dat";
            std::ofstream f_tmp;
            f_tmp.open(dataFileName.str().c_str(), std::ios::out);
            double* p_magnetization = magnetization;
            utilities::LatticeVector2D<double> kVector(offsetX,offsetY);
            utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1, (double)1/dimX);
            utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2, (double)1/dimY);
            f_tmp<<dimX<<"\n";
            f_tmp<<dimY<<"\n";
            f_tmp<<std::setw(15)<<std::left<<"kx values"<<"\t"<<std::setw(15)<<std::left<<"ky values"<<"\t"<<std::setw(15)<<std::left<<"magnetization"<<"\n";
            for(iSize_t x=0; x<dimX; ++x, kVector-=reciprocalLattice::G2, kVector += increment1)
            {
                for(iSize_t y=0; y<dimY; ++y, ++p_magnetization, kVector += increment2)
                {
                    f_tmp<<std::setw(15)<<std::left<<kVector.GetKx()<<"\t"<<std::setw(15)<<std::left<<kVector.GetKy()<<"\t"<<std::setw(15)<<std::left<<*(p_magnetization)<<"\n";
                }
            }
            f_tmp.close();
            utilities::cout.SecondaryOutput()<<"\n\tData output to file "<<dataFileName.str()<<std::endl;
            bool plotMagnetization;
            GetOption(optionList, plotMagnetization, "plot-magnetization", _LINE_, mpi);
            if(plotMagnetization)
            {
                utilities::cout.AdditionalInfo()<<"\n\tMAKING PLOT..."<<std::endl;
                std::string outPath;
                GetOption(optionList, outPath, "out-path", _LINE_, mpi);
                std::stringstream plotFileName;
                plotFileName.str("");
                plotFileName << outPath << "/magnetization_map_kx_" << dimX << "_ky_" << dimY << ".pdf";
                std::stringstream pythonScript; 
                pythonScript.str();
                pythonScript<<"#! //usr//bin//env python"<<_PYTHON_VERSION_<<"\n"
                "import matplotlib                              \n"
                "import numpy as np                             \n"
                "import math                                    \n"
                "matplotlib.use('Agg')                        \n\n"
                "import matplotlib.pyplot as plt              \n\n"
                "from mpl_toolkits.mplot3d import Axes3D      \n\n"
                "xData = []                                     \n"
                "yData = []                                     \n"   
                "magnetization = []                           \n\n"
                "fin=open(\""<<dataFileName.str()<<"\")       \n\n"
                "fin.readline()                                 \n"
                "fin.readline()                                 \n"
                "fin.readline()                                 \n"
                "for line in fin:                               \n"
                "    line = line.split()                        \n"
                "    xData.append(float(line[0]))               \n"
                "    yData.append(float(line[1]))               \n"
                "    magnetization.append(float(line[2]))       \n"
                "fin.close()                                  \n\n"
                "fig = plt.figure()                           \n\n"
                "ax = fig.gca(projection='3d')                \n\n"
                "ax.plot_trisurf(xData,yData,magnetization)     \n"
                "plt.savefig(\""<<plotFileName.str()<<"\") \n"
                "plt.close() \n";
                utilities::Script myScript;
                myScript.SetScript(pythonScript.str());
                myScript.Execute();
                utilities::cout.MainOutput()<<"\n\t...PLOT SUCESSFULLY GENERATED"<<std::endl;
                utilities::cout.MainOutput()<<"\n\tFind it at "<<plotFileName.str()<<std::endl;
            }
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate the set of coefficients and phase factors required to 
    //! change the single particle orbital basis to the basis of hybrid
    //! localized Wannier states as defined in Eq 32 of PRB 86,085129 (2012)
    //!
    //! The change of basis can be written as:
    //!  
    //! |x,k_y> = exp(-2*I*PI*offsetX x/dimX) exp^{Iphi_y(x,ky)}/sqrt(dimX) 
    //! sum_kx exp(-2*I*PI*kx x/dimX) {[[lambda_x(kx)]^kx]/product_{kappa=0}^{kx} 
    //! A_x (kappa,ky)} |kx,ky>
    //!
    //! The change of basis can be encoded in a matrix:
    //!
    //! c^+_{x,ky} = sum_{kx,ky} U_{(x,ky),(kx,ky)} c^+_{kx,ky}
    //!
    //! OR
    //!
    //! c^+_{kx,ky} = sum_{x,ky} U^-1_{(kx,ky),(x,ky)} c^+_{x,ky}
    //!
    //! where U is unitary and diagonal in ky. Then
    //!
    //! V_{(x,ky)_1,...,(x,ky)_4} = sum_{(kx,ky)_1,...,(kx,ky)_4)} V_{(kx,ky)_1,...,(kx,ky)_4}
    //!                           * U^-1_{(kx,ky)_1,(x,ky)_1} U^-1_{(kx,ky)_2,(x,ky)_2} 
    //!                           * U_{(x,ky)_3,(kx,ky)_3} U_{(x,ky)_4,(kx,ky)_4}
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModelGrid::GenerateWannierCoefficients(
        const iSize_t dimX,             //!<    X-dimension of k-space grid
        const iSize_t dimY,             //!<    Y-dimension of k-space grid
        const double offsetX,           //!<    Offset in the kx grid
        NonInteractingOflModel* blochTable,
                                        //!<    Table of Bloch wave functions
                                        //!     corresponding to that k-space grid
        dcmplx* M,                      //!<    A dimX*dimY by dimX*dimY matrix encoding the 
                                        //!     change of basis (to be populated)
        const utilities::MpiWrapper& mpi)   
                                        //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<"\n\t- CALCULATING CHANGE OF BASIS COEFFICIENTS";
            fflush(stdout);
        }
        //  First, we generate all of the data that will be needed to calculate the
        //  phase factor phi and the coefficients in the change of basis definition
        //////  Here we implement Eq 26     ////////////////////////////////////////////
        dcmplx* Ax = new dcmplx[dimX*dimY];
        for(kState_t kx=0; kx<dimX; ++kx)
        {
            for(kState_t ky=0; ky<dimY; ++ky)
            {
                kState_t kxPlusOnePeriodic = kx==(dimX-1) ? 0 : kx+1;
                dcmplx aTempX = utilities::linearAlgebra::DotProduct<dcmplx>(blochTable[kx*dimY+ky].m_blochCoefficients,
                                blochTable[kxPlusOnePeriodic*dimY+ky].m_blochCoefficients, blochTable->m_dim);
                //  Set the Berry connection to the normalized accumulated value
                //  as defined in Eq. 27:
                Ax[kx*dimY+ky] = aTempX/abs(aTempX);
            }
        }
        //  As defined in Eq. 30:
        dcmplx* lambdaX = new dcmplx[dimY];
        for(kState_t ky=0; ky<dimY; ++ky)
        {
            lambdaX[ky] = 1.0;
        }
        //  Generate products of Ax values
        for(kState_t kx=0; kx<dimX; ++kx)
        {
            for(kState_t ky=0; ky<dimY; ++ky)
            {
                lambdaX[ky] *= Ax[kx*dimY+ky];
            }
        }
        for(kState_t ky=0; ky<dimY; ++ky)
        {
            lambdaX[ky] = std::pow(lambdaX[ky], 1.0/dcmplx(dimX,0));
        }
        //  aProd(kx,ky) = product_{kappa=0}^{kx-1} Ax(kappa,ky)
        dcmplx* aProd = new dcmplx[dimX*dimY];
        for(kState_t ky=0; ky<dimY; ++ky)
        {
            aProd[ky] = 1.0;
        }
        for(kState_t kx=1; kx<dimX; ++kx)
        {
            for(kState_t ky=0; ky<dimY; ++ky)
            {
                aProd[kx*dimY+ky] = Ax[(kx-1)*dimY+ky]*aProd[(kx-1)*dimY+ky];
            }
        }
        delete[] Ax;    //<-- move these delete[] commands if implementing below!!
        //delete[] Ay;
        
        ////////////////////////////////////////////////////////////////////////////////
        #if _NOT_IMPLEMENTED_
        
        //  THE GAUGE-FIXED PART OF THE WANNIER FUNCTION PERSCRIPTION IS NOT CURRENTLY IMPLEMETED
        dcmplx  lambdaY = 1.0;
        for(kState_t ky=0;ky<dimY;++ky)
        {
            lambdaY *= Ay[ky];
        }
        lambdaY = std::pow(lambdaY, 1.0/(double)dimY);
        //////  Tabulate U(ky) as defined in Eq 53      ////////////////////////////////
        dcmplx* U = new dcmplx[dimY];
        for(kState_t ky=0; ky<dimY; ++ky)
        {
            U[ky] = 0.0;
        }
        for(kState_t kx=0; kx<dimX; ++kx)
        {
            for(kState_t ky=0; ky<dimY; ++ky)
            {
                kState_t kyPlusOnePeriodic = (ky==(dimY-1) ? 0 : ky+1);
                dcmplx W = Ay[kx*dimY+ky]*aProd[kx*dimY+ky]/(aProd[kx*dimY+kyPlusOnePeriodic]*Ay[ky]);
                dcmplx Wbar = std::pow((lambdaX[ky]/lambdaX[kyPlusOnePeriodic]), (double)kx);
                U[ky] += W/Wbar;
            }
        }
        for(kState_t ky=0; ky<dimY; ++ky)
        {
            U[ky] /= abs(U[ky]);
        }
        //////      Calculate and store omega, defined in Eq 59     ////////////////////
        //  omega = [product_ky U(ky)]^(1/dimY)  
        dcmplx omega = 1.0;
        for(kState_t ky=0;ky<dimY;++ky)
        {
            omega *= U[ky];
        }
        omega = std::pow(omega, 1.0/(double)dimY);
        //////  Fix the phase factors phi_y(x,ky)       ////////////////////////////////
        //  By definition, the 0,0 phase factor is set to 0 (so exp(phase) = 1.0). 
        //  The remaining factors are calculated from Eq 63.
        dcmplx* phaseFactors = new dcmplx[dimY];
        phaseFactors[0] = 1.0;
        //  Work out the remaining pre-factors from that one
        for(kState_t ky=0; ky<dimY-1; ++ky)
        {
            phaseFactors[ky+1] = phaseFactors[ky]*lambdaY*omega/(Ay[ky]*U[ky]);
        }
        delete[] U;             //<-- move these delete[] commands if implemented!!
        delete[] phaseFactors;
        #endif
        ////////////////////////////////////////////////////////////////////////////////
        
        //////      Now populate the change of basis matrix     ////////////////////
        iSize_t dim = dimX*dimY;
        for(kState_t x=0; x<dimX; ++x)
        {
            for(kState_t ky=0; ky<dimY; ++ky)
            {
                iSize_t index1 = x*dimY+ky;
                for(kState_t kx=0; kx<dimX; ++kx)
                {
                    for(kState_t ky1=0; ky1<dimY; ++ky1)
                    {
                        iSize_t index2 = kx*dimY+ky1;
                        if(ky == ky1)
                        {
                            M[index1*dim+index2] = exp(-2.0*I*PI*((double)kx+offsetX)*(double)x/(double)dimX)*std::pow(lambdaX[ky], dcmplx(kx,0))/(sqrt(dimX)*aProd[index2]);
                        }
                        else
                        {
                            M[index1*dim+index2] = 0.0;
                        }
                    }
                }
            }
        }    
        //////      Clear memory that we no longer need     ////////////////////////////
        delete[] lambdaX;
        delete[] aProd;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<" - DONE"<<std::endl;
        }  
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to obtain a map of the sigma^z values associated with each
    //!  k-space point
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModelGrid::CalcualteSigmaZMap(
        const iSize_t dimX,             //!<    X-dimension of k-space grid
        const iSize_t dimY,             //!<    Y-dimension of k-space grid
        NonInteractingOflModel* blochTable,
                                        //!<    Table of Bloch wave functions
                                        //!     corresponding to that k-space grid
        double* map)                    //!<    A dimX*dimY array to store the sigma^z_k 
                                        //!     values
    {
        double* p_map = map;
        for(kState_t kx=0; kx<dimX; ++kx)
        {
            for(kState_t ky=0; ky<dimY; ++ky, ++p_map)
            {
                *(p_map) = blochTable[kx*dimY+ky].ApplySigmaZ(0);
            }
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to generate the single-particle orbitals
    //!
    //! psi^sigma_{kx,ky} (r) = |<r,sigma|kx,ky>|^2
    ////////////////////////////////////////////////////////////////////////////////
    void NonInteractingOflModelGrid::CalculateSpatialWaveFunctions(
        boost::program_options::variables_map* optionList,
                                            //!<    Command line arguments
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t CALCULATING MAPPING TO SPATIAL WAVE FUNCTIONS"<<std::endl;
        }
        NonInteractingOflModelData params = NonInteractingOflModelData(optionList, mpi);          
        double matrixElementTol;
        iSize_t minBlochCutOff;
        bool useWannierBasis;
        iSize_t dimX;
        iSize_t dimY;
        iSize_t realGridDim;
        double offsetX;
        double offsetY;
        double gridSpacing;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
	        int getBasis;
	        GetOption(optionList, matrixElementTol, "tol", _LINE_, mpi);
	        GetOption(optionList, realGridDim, "calculate-spatial-wavefunctions", _LINE_, mpi);
	        GetOption(optionList, dimX, "x-grid", _LINE_, mpi);
	        GetOption(optionList, dimY, "y-grid", _LINE_, mpi);
	        GetOption(optionList, getBasis, "sptial-basis", _LINE_, mpi);
	        GetOption(optionList, minBlochCutOff, "x-cut", _LINE_, mpi);
	        GetOption(optionList, offsetX, "x-grid-shift", _LINE_, mpi);
	        GetOption(optionList, offsetY, "y-grid-shift", _LINE_, mpi);
	        GetOption(optionList, gridSpacing, "grid-spacing", _LINE_, mpi);
            if(1 == getBasis)
            {
                useWannierBasis = true;
            }
            else
            {
                useWannierBasis = false;
            }
        }
        mpi.ExitFlagTest();
        mpi.Sync(&matrixElementTol, 1, 0);
        mpi.Sync(&realGridDim, 1, 0);
        mpi.Sync(&useWannierBasis, 1, 0);
        mpi.Sync(&dimX, 1, 0);
        mpi.Sync(&dimY, 1, 0);
        mpi.Sync(&minBlochCutOff, 1, 0);
        mpi.Sync(&offsetX, 1, 0);
        mpi.Sync(&offsetY, 1, 0);
        mpi.Sync(&gridSpacing, 1, 0);
        params.MpiSynchronize(0, mpi);
        if((int)(realGridDim*realGridDim)<mpi.m_nbrProcs)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
	        {
	            std::cerr<<"\n\tERROR IN CalculateSpatialWaveFunctions: number of nodes exceeds grid size"<<std::endl;
	        }
	        exit(EXIT_FAILURE);
        }
        //////      First generate the Bloch wave functions     ////////////////////
        NonInteractingOflModel* blochWaveFunctionTable = new NonInteractingOflModel[dimX*dimY]; 
        this->GenerateBlochWaveFunctionTable(minBlochCutOff, dimX, dimY, offsetX, offsetY, blochWaveFunctionTable, &params, mpi);
        //////      If required, calculate the change of basis matrix       ////////
        dcmplx* basisChangeMatrix = 0;
        if(useWannierBasis)
        {
            //  Allocate space to store the change of basis matrix U_{(x,ky),(kx,ky)}
            basisChangeMatrix = new dcmplx[dimX*dimY*dimX*dimY];
            this->GenerateWannierCoefficients(dimX, dimY, offsetX, blochWaveFunctionTable, basisChangeMatrix, mpi);
        }
        //////      Convert to spatial wave functions on a real space grid     /////
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<"\n\t- CONVERTING TO SPATIAL WAVE FUNCTION GRID"<<std::endl;
        }
        //  A container to store the map between the Bloch function basis and the associated
        //  real-space wave functions psi^sigma_{kx,ky} (r)
        utilities::MultiHashMap<double> spatialWaveFunctionTable;
        //  For parallel processes, divide up tasks amongst available nodes
        mpi.DivideTasks(mpi.m_id, realGridDim*realGridDim, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
        iSize_t counter  = 0;
        utilities::LoadBar progress;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            progress.Initialize(mpi.m_lastTask - mpi.m_firstTask+1);
        }
        for(kState_t rx=0; rx<realGridDim; ++rx)
        {
            for(kState_t ry=0; ry<realGridDim; ++ry, ++counter)
            {  
                if(mpi.m_firstTask <= (int)counter && mpi.m_lastTask >= (int)counter)
                {
                    if(useWannierBasis)
                    {
                        //  If we are using the Wannier basis then we need to include the change
                        //  of basis in order to plot the spatial wave functions. The wave 
                        //  functions are given by
                        //
                        //  psi^sigma_{x,ky} (r) = exp(-2*I*PI*offsetX x/N_x) sum_{k_x=0}^{N_x-1} exp(-2*I*pi*k_x*x/N_x) 
                        //  {[[lambda_x(kx)]^kx]/product_{kappa=0}^{kx} A_x (kappa,ky)} sum_G a^{kx,ky}_{G,sigma} exp^(ir.(G+{kx,ky})
                        //  This can be achieved by performing a set of N_y 1D Fourier transforms
                        //  Allocate memory to store input/outputs of Fourier transforms
                        dcmplx* fourierCoefficientsUp   = new dcmplx[dimX*dimY];
                        dcmplx* outputCoefficientsUp    = new dcmplx[dimX*dimY];
                        dcmplx* fourierCoefficientsDown = new dcmplx[dimX*dimY];
                        dcmplx* outputCoefficientsDown  = new dcmplx[dimX*dimY];
                        //  Populate the Fourier transform inputs
                        NonInteractingOflModel*  p_blochTable = blochWaveFunctionTable;
                        for(kState_t ky=0; ky<dimY; ++ky)
                        {
                            for(kState_t kx=0; kx<dimX; ++kx, ++p_blochTable)
                            {
                                //  This value is given by
                                //  {[[lambda_x(kx)]^kx]/product_{kappa=0}^{kx} A_x (kappa,ky)} sum_G a^{kx,ky}_{G,sigma} exp^(ir.(G+{kx,ky})
                                //  Note that the values stored in the change of basis matrix are
                                //  U[x,ky,kx,ky] = exp(-2*I*PI*offsetX x/N_x) exp(-2*I*pi*k_x*x/N_x)/sqrt(N_x)*{[[lambda_x(kx)]^kx]/product_{kappa=0}^{kx} A_x (kappa,ky)}
                                //  Since we don't need the complex exponential part,
                                //  (as it will be implemted by the Fourier transform)
                                //  we need only extract U[0,ky,kx,ky]*sqrt(N_x)
                                
                                unsigned int index = (ky)*dimX*dimY+(ky*dimX+kx);
                                fourierCoefficientsUp[ky*dimX+kx]   = basisChangeMatrix[index]*sqrt(dimX)*p_blochTable->GetSpatialWaveFunction(_SPIN_UP_, rx, ry, gridSpacing);
                                fourierCoefficientsDown[ky*dimX+kx] = basisChangeMatrix[index]*sqrt(dimX)*p_blochTable->GetSpatialWaveFunction(_SPIN_DOWN_, rx, ry, gridSpacing);
                            }
                        }
                        //  Perform the set of transforms
                        utilities::DiscreteFourierTransform1D(dimX, dimY, fourierCoefficientsUp, outputCoefficientsUp, FFTW_BACKWARD);
                        utilities::DiscreteFourierTransform1D(dimX, dimY, fourierCoefficientsDown, outputCoefficientsDown, FFTW_BACKWARD);
                        //  Deallocate memory for Fourier transforms
                        delete[] fourierCoefficientsUp;
                        delete[] fourierCoefficientsDown;
                        //  Store the set of coefficients psi^sigma_{x,ky} (r)
                        dcmplx* p_outputCoefficientsUp = outputCoefficientsUp;
                        dcmplx* p_outputCoefficientsDown = outputCoefficientsDown;
                        kState_t combinedR = this->CombineIndexes(rx, ry, realGridDim);
                        for(kState_t ky=0; ky<dimY; ++ky)
                        {
                            for(kState_t x=0; x<dimX; ++x, ++p_outputCoefficientsUp, ++p_outputCoefficientsDown)
                            {         
                                kState_t combinedK = this->CombineIndexes(x, ky, dimY);
                                spatialWaveFunctionTable.Insert(utilities::Key(combinedK, combinedR, (bool)_SPIN_UP_)) = abs(*p_outputCoefficientsUp)*abs(*p_outputCoefficientsUp);
                                spatialWaveFunctionTable.Insert(utilities::Key(combinedK, combinedR, (bool)_SPIN_DOWN_)) = abs(*p_outputCoefficientsDown)*abs(*p_outputCoefficientsDown);
                            }
                        }
                        //  Deallocate memory for Fourier transform outputs
                        delete[] outputCoefficientsUp;
                        delete[] outputCoefficientsDown;
                    }
                    else
                    {
                        //  If we use the regular Bloch wave function basis, then the spatial
                        //  wave functions are given by the sums of Bloch coefficients we've
                        //  already calculated
                        //  psi^sigma_{kx,ky} (r) = sum_G a^{kx,ky}_{G,sigma} exp^(ir.(G+{kx,ky})
                        kState_t combinedR = this->CombineIndexes(rx, ry, realGridDim);
                        NonInteractingOflModel*  p_blochTable = blochWaveFunctionTable;
                        for(kState_t ky=0; ky<dimY; ++ky)
                        {
                            for(kState_t kx=0; kx<dimX; ++kx, ++p_blochTable)
                            {            
                                kState_t combinedK = this->CombineIndexes(kx, ky, dimY);
                                double absValueUp = abs(p_blochTable->GetSpatialWaveFunction(_SPIN_UP_, rx, ry, gridSpacing));
                                double absValueDown = abs(p_blochTable->GetSpatialWaveFunction(_SPIN_DOWN_, rx, ry, gridSpacing));
                                spatialWaveFunctionTable.Insert(utilities::Key(combinedK, combinedR, (bool)_SPIN_UP_)) = absValueUp;
                                spatialWaveFunctionTable.Insert(utilities::Key(combinedK, combinedR, (bool)_SPIN_DOWN_)) = absValueDown;
                            }
                        }
                    }
                    if(0 == mpi.m_id && mpi.m_lastTask >= (int)counter)	// FOR THE MASTER NODE
                    {
                        progress.Display(counter+1);
                    }
                }
            }
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<std::endl;
        }
        //////      Output spatial wave function data to a text file        ////////        
        std::stringstream dataFileName;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {   
	        std::string outPath;
	        GetOption(optionList, outPath, "out-path", _LINE_, mpi);
            dataFileName.str("");
            dataFileName << outPath << "/spatial_wave_function_kx_" << dimX << "_ky_" << dimY;
            bool useSql;
            GetOption(optionList, useSql, "use-sql", _LINE_, mpi);
            if(useSql)
            {
                iSize_t sqlId;
                GetOption(optionList, sqlId, "sql-id", _LINE_, mpi);
                dataFileName << "_id_" << sqlId;
            }
            dataFileName << "_grid_size_" << realGridDim;
            if(useWannierBasis)
            {
                dataFileName << "_wannier_basis";
            }
            dataFileName << ".dat";
        }
        io::fileFormat_t format = io::_TEXT_;
        std::ofstream f_out;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            GenFileStream(f_out, dataFileName.str(), format, mpi);
        }
        mpi.ExitFlagTest();
        spatialWaveFunctionTable.ToFile(f_out, format, 3, mpi);
        if(f_out.is_open())
        {
            f_out.close();
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {   
            utilities::cout.SecondaryOutput()<<"\n\tFinished writing spatial wave funciton values to a file. See the file here: "<<dataFileName.str()<<std::endl;
        }
        return;
    }
}   //  End namespace diagonalization
