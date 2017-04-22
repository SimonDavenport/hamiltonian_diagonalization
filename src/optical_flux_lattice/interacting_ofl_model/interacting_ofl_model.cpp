////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!	 \file
//!     This file defines a class to store the optical flux lattice interacting
//!     Hamiltonian. See e.g. PRL 109, 265301 (2012)	           
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
#include "interacting_ofl_model.hpp"

namespace diagonalization
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Default constructor for the InteractingOflModel class
    //! 
    //! NOTE: implicitly calls the default constructors of all internal classes
    //!
    ////////////////////////////////////////////////////////////////////////////////
    InteractingOflModel::InteractingOflModel()
    {}

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor for the InteractingOflModel class
    //!
    //! Sets the internal instance of the hamiltonianData to an external value.
    //!
    ////////////////////////////////////////////////////////////////////////////////
    InteractingOflModel::InteractingOflModel(
        InteractingOflModelData data) //!<    External instance of InteractingOflModelData
        :   m_params(data)
    {}

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Constructor from command line arguments 
    ////////////////////////////////////////////////////////////////////////////////
    InteractingOflModel::InteractingOflModel(
        boost::program_options::variables_map* optionList,    
                                                //!<    Parsed command line argument list
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        //  Set model parameters based on command line arguments
        m_params = InteractingOflModelData(optionList, mpi);
        m_oflParameters = NonInteractingOflModelData(optionList, mpi);
	    if(m_params.m_useWannierBasis)
        {
            //  Reserve memory to store basis change matrix
            m_basisChangeMatrix.resize(m_params.m_dimX*m_params.m_dimY*m_params.m_dimX*m_params.m_dimY);
        }
        //  Set ARPACK related options
        m_hamiltonian.InitializeArpack(optionList, mpi);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Destructor for the InteractingOflModel class
    ////////////////////////////////////////////////////////////////////////////////
    InteractingOflModel::~InteractingOflModel()
    {}

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Update the status of the SQL database once calculations have
    //! been completed
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::UpdateSqlStatus(
        const utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        //  Update the SQL status to record any completed operations
        m_params.UpdateSqlStatus(mpi);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Set the kxTot and kyTot linear momentum values. If this function is 
    //! called then the class automatically switches to block diagonalization mode,
    //! otherwise all sectors are diagonalized.
    //!
    //! If the Wannier basis is used, then the kx linear momentum quantum number
    //! will be ignored. 
    ////////////////////////////////////////////////////////////////////////////////

    void InteractingOflModel::SetSector(
        const kState_t kxSector,     //!<    Specify kx total linear momentum sector
        const kState_t kySector)     //!<    Specify ky total linear momentum sector
    {
        m_params.m_kxTot = kxSector;
        m_params.m_kyTot = kySector;
        m_params.m_blockDiagonalize = true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate tables describing the quadratic and quartic terms in the
    //! Hamiltonian
    //!
    //! This implements the Hamiltonian given e.g. in eq (5) of PRL 109,265301
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::BuildTermTables(
        boost::program_options::variables_map* optionList,
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        //////  POPULATE LOOK-UP TABLES FOR THE QUADRATIC AND QUARTIC TERM  //////////
        /////   COEFFICIENTS OF OUR HAMILTONIAN                             //////////
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    { 
            utilities::cout.MainOutput()<<"\n\t============ BUILDING MATRIX ELEMENT LOOK-UP TABLES ============"<<std::endl;
        }
        double matrixElementTol;
        iSize_t minBlochCutOff;
        iSize_t maxBlochCutOff;
        iSize_t cutOffStep;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            GetOption(optionList, matrixElementTol, "tol-vkkkk", _LINE_, mpi);
            GetOption(optionList, minBlochCutOff, "min-k-cut", _LINE_, mpi);
            GetOption(optionList, maxBlochCutOff, "max-k-cut", _LINE_, mpi);
            GetOption(optionList, cutOffStep, "k-cut-step", _LINE_, mpi);
        }
        mpi.ExitFlagTest();
        mpi.Sync(&matrixElementTol, 1, 0);
        mpi.Sync(&minBlochCutOff, 1, 0);
        mpi.Sync(&maxBlochCutOff, 1, 0);
        mpi.Sync(&cutOffStep, 1, 0);
        //  Calculate the Vkkkk coefficients iteratively, increasing the value of the
        //  cut off in the Bloch coefficient calculation until we converge on a
        //  stable result to within some tolerance level
        //  we also use these solutions to populate the quadratic term table
        m_quarticTables.Initialize(m_params.m_dimX*m_params.m_dimY, mpi);
        std::vector<int> gxTable(m_quarticTables.GetDimension());
        std::vector<int> gyTable(m_quarticTables.GetDimension());
        this->Generate4KTable(m_quarticTables.GetKTable(), &gxTable, &gyTable, mpi);
        //  Allocate memory to store a temporary look-up table of quartic terms
        //  but only store on the master node
        std::vector<dcmplx> listOfQuarticTerms;
        if(0 == mpi.m_id) // FOR THE MASTER NODE
        {   
            listOfQuarticTerms.resize(m_quarticTables.GetDimension());
        }
        //    Allocate memory to store a list of Bloch wave function tables
        NonInteractingOflModel* blochWaveFunctionTable = new NonInteractingOflModel[m_params.m_dimX*m_params.m_dimY]; 
        iSize_t cutOff = minBlochCutOff;
        double currentTol = 1000.0;
        double prevTol = currentTol;
        //  Generate the first approximation
        NonInteractingOflModelGrid modelGrid;
        modelGrid.GenerateBlochWaveFunctionTable(cutOff, m_params.m_dimX, m_params.m_dimY, m_params.m_offsetX, 
                                                  m_params.m_offsetY, blochWaveFunctionTable, &m_oflParameters, mpi);
        this->GenerateQuarticTerms(m_quarticTables.GetVTable(), m_quarticTables.GetDimension(), blochWaveFunctionTable, 
                                   &gxTable, &gyTable, mpi);
        do
        {
            if(0 == mpi.m_id) // FOR THE MASTER NODE
            {
                if(cutOff>maxBlochCutOff)
                {
                    std::cerr<<"WARNING in BuildTermTables: max Bloch cut off reached without converging!"<<std::endl;
                    break;
                }
                utilities::cout.DebuggingInfo()<<"\n\t- PERFORMING MORE ACCURATE CALCULATION "<<std::endl;
            }
            mpi.ExitFlagTest();
            //  Increase the range of Bloch coefficients used in the calculation
            cutOff += cutOffStep;
            //  Re-generate the coefficient table
            modelGrid.GenerateBlochWaveFunctionTable(cutOff, m_params.m_dimX, m_params.m_dimY, m_params.m_offsetX, 
                                                      m_params.m_offsetY, blochWaveFunctionTable, &m_oflParameters, mpi);
            this->GenerateQuarticTerms(&listOfQuarticTerms, m_quarticTables.GetDimension(), blochWaveFunctionTable, 
                                       &gxTable, &gyTable, mpi);
            //  Get the maximum abs difference in Vkkkk values 
            if(0 == mpi.m_id) // FOR THE MASTER NODE
            {
                auto it_list1 = listOfQuarticTerms.begin();
                auto it_list2 = m_quarticTables.GetVTable()->begin();
                currentTol = 0.0;
                for(iSize_t i=0; i<m_quarticTables.GetDimension(); ++i, ++it_list1, ++it_list2)
                {
                    const double temp = fabs(abs(*it_list1)-abs(*it_list2));
                    if(temp>currentTol)
                    {
                        currentTol = temp;
                    }
                }
                *(m_quarticTables.GetVTable()) = listOfQuarticTerms;
                utilities::cout.SecondaryOutput()<<"\n\tBLOCH BAND CUT OFF: "<<cutOff<<"\tMAX CHANGE IN |Vkkkk|: "<<currentTol<<std::endl;
            }
            mpi.Sync(&currentTol, 1, 0);
            //  Break out if the previous tolerance was better than the current 
            //  tolerance (a precision limit has been reached)
            if(fabs(prevTol)<fabs(currentTol))
            {
                std::cerr<<"\n\tWARNING: DUE TO NON CONVERGENCE, Vkkkk CALCULATION EXITED BEFORE SPECIFIED TOLERENCE REACHED!"<<std::endl;
                break;
            }
            else
            {
                prevTol = currentTol;
            }
        }
        while(currentTol>matrixElementTol);
        //  Generate the corresponding table of quadratic terms
        m_quadraticTables.Initialize(m_params.m_dimX*m_params.m_dimY,mpi);
        if(0 == mpi.m_id)    // FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<"\n\t- MAXIMUM TOLERENCE ("<<std::max(matrixElementTol,fabs(currentTol))<<") REACHED!"<<std::endl;
            auto it_quadraticTable = m_quadraticTables.GetVTable()->begin();
            for(iSize_t i=0; i<m_quadraticTables.GetDimension(); ++i, ++it_quadraticTable)
            { 
                *(it_quadraticTable) = blochWaveFunctionTable[i].GetBandEnergy(0);
            }
        }
        //  Multiply all quartic table values by a constant coefficient
        if(0 == mpi.m_id)    // FOR THE MASTER NODE
        {   
            m_quarticTables.SetCoefficient(0.5*m_params.m_interactionStrength);
        }
        //  Synchronize with the master node
        if(0 == mpi.m_id)    // FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<"\n\t- MPI SYNC HAMILTONIAN COEFFICIENT TABLES"<<std::endl;
        }
        m_quadraticTables.MpiSynchronize(0, mpi);
        m_quarticTables.MpiSynchronize(0, mpi);          
        //  Calculate sigmaZ map
        this->CalcualteSigmaZMap(blochWaveFunctionTable, mpi);
        //  Convert to Wannier basis if requested
        if(m_params.m_useWannierBasis)
        {
            this->ChangeToWannierBasis(blochWaveFunctionTable, matrixElementTol, mpi);
        }
        //  g Tables and blochWaveFunctionTable no longer required
        delete[] blochWaveFunctionTable;
        if(_HASH_ == m_params.m_setTableFormat)
        {
            this->ConvertTableFormat(mpi);
        }
        m_params.m_termTablesBuilt = true;
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the unique k1 value corresponding to a given k2,k3
    //! and k4 assuming 2D momentum conservation.
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::Generate4KTable(
        std::vector<kState_t>* kTable,      //!<    Preallocated table of k-values
        std::vector<int>* gxTable,          //!<    Preallocated table of reciprocal lattice values
        std::vector<int>* gyTable,          //!<    Preallocated table of reciprocal lattice values
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING MOMENTUM CONSERVING k1,k2,k3,k4 TABLE"<<std::endl;
        }
        kState_t kMax = m_params.m_dimX*m_params.m_dimY;
        auto it_kTable = kTable->begin();
        auto it_gxTable = gxTable->begin();
        auto it_gyTable = gyTable->begin();
        for(kState_t k4=0; k4<kMax; ++k4)
        {
            for(kState_t k3=k4; k3<kMax; ++k3)
            {
                for(kState_t k2=0; k2<kMax; ++k2, ++it_kTable, ++it_gxTable, ++it_gyTable)
                {
                    if(it_kTable < kTable->end() && it_gxTable < gxTable->end() && it_gyTable < gyTable->end())
                    {
                        //  For every k2,k3,k4 there is a unique k1 in the BZ
                        //  such that k1+k2-k3-k4 = a reciprocal lattice vector 
                        //  as G = x GX + y GY
                        int x,y;
                        kState_t k1 = this->FindK1(k2, k3, k4, x, y);
                        *(it_kTable) = k1;
                        *(it_gxTable) = x;
                        *(it_gyTable) = y;
                    }
                }
            }
        }
	    MPI_Barrier(mpi.m_comm);
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the unique k1 value corresponding to a given k2,k3
    //! and k4 assuming 2D momentum conservation (Wannier basis)
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::GenerateWannier4KTable(
        utilities::MultiHashMultiMap<kState_t>* kHashTable,  
                                            //!<    Hash table container
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING MOMENTUM CONSERVING k1,k2,k3,k4 TABLE"<<std::endl;
        }
        for(kState_t kx4=0; kx4<m_params.m_dimX; ++kx4)
        {
            for(kState_t ky4=0; ky4<m_params.m_dimY; ++ky4)
            {
                for(kState_t kx3=0; kx3<m_params.m_dimX; ++kx3)
                {
                    for(kState_t ky3=0; ky3<m_params.m_dimY; ++ky3)
                    {
                        for(kState_t kx2=0; kx2<m_params.m_dimX; ++kx2)
                        {
                            for(kState_t ky2=0; ky2<m_params.m_dimY; ++ky2)
                            {
                                //  Determine an appropriate ky1 value
                                //  given linear momentum conservation
                                //  along the ky direction
                                int gTot2;
                                kState_t ky1 = this->FindK1(ky2, ky3, ky4, gTot2);
                                for(kState_t kx1=0; kx1<m_params.m_dimX; ++kx1)
                                {
                                    kHashTable->Insert(kx1*m_params.m_dimY+ky1, utilities::Key(kx2*m_params.m_dimY+ky2, 
                                                       kx3*m_params.m_dimY+ky3, kx4*m_params.m_dimY+ky4));
                                }
                            }
                        }
                    }
                }
            }
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the unique k1 value corresponding to a given k2
    //! assuming Wannier basis is used
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::GenerateWannier2KTable(
        utilities::MultiHashMultiMap<kState_t>* kHashTable,   //!<    Hash table container
        const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING MOMENTUM CONSERVING k1,k2 TABLE"<<std::endl;
        }
        //  Only ky linear momentum is diagonal
        for(kState_t kx1=0; kx1<m_params.m_dimX; ++kx1)
        {
            for(kState_t ky1=0; ky1<m_params.m_dimY; ++ky1)
            {
                for(kState_t kx2=0; kx2<m_params.m_dimX; ++kx2)
                {
                    kHashTable->Insert(kx1*m_params.m_dimY+ky1, kx2*m_params.m_dimY+ky1);
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Determine the value k1 such that k3+k4-k1-k2 + G = 0 with k1,k2,k3,k4
    //!  in the BZ and G a reciprocal lattice vector of the hexagonal lattice.
    //!
    //! \return an integer index defining the location of k1 within a 2D BZ grid
    ////////////////////////////////////////////////////////////////////////////////
    kState_t InteractingOflModel::FindK1(
        const kState_t k2,        //!<   Index of k2 vector in the BZ
        const kState_t k3,        //!<   Index of k3 vector in the BZ
        const kState_t k4,        //!<   Index of k4 vector in the BZ
        int& gTot1,               //!<   Counter for number G1 lattice vectors shifted
        int& gTot2)               //!<   Counter for number G2 lattice vectors shifted
        const
    {
        //  Starting value at k1 = k3 + k4 - k2;
        //  Convert k3,k2,k1 indices to a 2D grid
        int finalKx = (int)floor((double)k3/(m_params.m_dimY)) + (int)floor((double)k4/(m_params.m_dimY)) - (int)floor((double)k2/(m_params.m_dimY));
        int finalKy = (int)utilities::Modulo<int>(k3, m_params.m_dimY) + (int)utilities::Modulo<int>(k4, m_params.m_dimY) - (int)utilities::Modulo<int>(k2, m_params.m_dimY);
        //  Determine the vector value shifted to the BZ
        int shiftedKx = utilities::Modulo<int>(finalKx, m_params.m_dimX);
        int shiftedKy = utilities::Modulo<int>(finalKy, m_params.m_dimY);
        // Determine the number of reciprocal lattice vectors shifted
        int gVectorX = shiftedKx - finalKx;
        int gVectorY = shiftedKy - finalKy;
        if(gVectorX < 0)
        {
            gTot1 = -floor((double) -gVectorX/m_params.m_dimX);
        }
        else
        {
            gTot1 = floor((double) gVectorX/m_params.m_dimX);
        }
        
        if(gVectorY < 0)
        {
            gTot2 = -floor((double) -gVectorY/m_params.m_dimY);
        }
        else
        {
            gTot2 = floor((double) gVectorY/m_params.m_dimY);
        }
        return shiftedKx*m_params.m_dimY+shiftedKy;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Determine a list of k1y values such that k3y+k4y-k1y-k2y + Gy = 0 
    //!
    //! \return The k1y index satisfying momentum conservation along the y direction
    ////////////////////////////////////////////////////////////////////////////////
    kState_t InteractingOflModel::FindK1(
        const kState_t k2y,        //!<   Index of k2y in the Wannier basis
        const kState_t k3y,        //!<   Index of k3y in the Wannier basis
        const kState_t k4y,        //!<   Index of k4y in the Wannier basis
        int& gTot2)                //!<   Counter for number G2 lattice vectors shifted
        const
    {
        //  Starting value at k1y = k3y + k4y - k2y;
        int finalKy = (int)k3y + (int)k4y - (int)k2y;
        //  Determine the value shifted to the BZ
        int shiftedKy = utilities::Modulo<int>(finalKy, m_params.m_dimY);
        // Determine the number of reciprocal lattice vectors shifted
        int gVector = shiftedKy - finalKy;
        if(gVector < 0)
        {
            gTot2 = -floor((double) -gVector/m_params.m_dimY);
        }
        else
        {
            gTot2 = floor((double) gVector/m_params.m_dimY);
        }
        return shiftedKy;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate the coefficients V_{k1,k2,k3,k4} in the interacting Hamiltonian
    //! 
    //! The coefficients are given by:
    //!
    //! \[ V_{k1,k2,k3,k4}  = \sum_\sigma \sum_{G_A,G_B,G_C,G_D} \delta (k3+k4-k1-k2+G_C+G_D-G_A-G_B)
    //! * c^{k1,*}_{G_A,\sigma} c^{k2,*}_{G_B,\bar{\sigma}} c^k3_{G_C,\bar{\sigma}} c^k4_{G_D,\sigma} \]
    //!
    //! The sum over G_A,...G_D is greatly simplified by writing is as two 
    //! independent sets of convolutions for the spin-up/spin-down terms separately.
    //!
    //! The calculation is optimized by mapping V_{k2,k1,k4,k3} to V_{k1,k2,k3,k4}
    //!
    //! In pseudo-code:
    //!
    //! V_up = {x*,y*,z*}^{x,y,z} = {xz*,xy*+yz*,xx*+yy*+zz*,yx*+zy*,zx*}
    //!
    //! The the sum over G_A,..G_D is given by V_up.V_down.
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::GenerateQuarticTerms(
        std::vector<dcmplx>* listOfQuarticTerms,     
                                        //!<    A pre-declared array to contain the
                                        //!     look-up table
        const iSize_t dimension,        //!<    Dimension of the quartic term table
                                        //!     (full value must be used on all nodes)
        NonInteractingOflModel* blochTable,
                                        //!<    Lookup table of Bloch wave functions
        std::vector<int>* gxTable,      //!<    Momentum conserving reciprocal lattice vectors (x)
        std::vector<int>* gyTable,      //!<    Momentum conserving reciprocal lattice vectors (y)
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
        const
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- GENERATING QUARTIC TERMS\n"<<std::endl;
        }
        //  Determine some dimensions used in calculating the convolutions
        const iSize_t blochTableXDim = 2*blochTable->GetCutOffX()+1;
        const iSize_t blochTableYDim = 2*blochTable->GetCutOffY()+1;
        const iSize_t convolutionDim = utilities::GetArrayConvolution2DSize(blochTableXDim, blochTableYDim);
        const iSize_t convolutionArrayOffset = utilities::GetArrayConvolution1DSize(blochTableXDim, blochTableYDim);
        //  Get the total number of k-states in the discretized Brillouin zone
        kState_t kMax = m_params.m_dimX*m_params.m_dimY;
        //  Allocate further memory to store temporary copies of the Bloch wave function tables
        //  and their convolutions
        dcmplx* blochAup     = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochAdown   = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochBup     = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochBdown   = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochCup     = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochCdown   = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochDup     = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* blochDdown   = new dcmplx[blochTableXDim*blochTableYDim];
        dcmplx* convolution1 = new dcmplx[convolutionDim];
        dcmplx* convolution2 = new dcmplx[convolutionDim];
        dcmplx* convolution3 = new dcmplx[convolutionDim];
        dcmplx* convolution4 = new dcmplx[convolutionDim];
        //  Populate the list of quadratic terms
        mpi.DivideTasks(mpi.m_id, dimension, mpi.m_nbrProcs, &mpi.m_firstTask, &mpi.m_lastTask, false);
        //  Allocate memory to store part of the result on each node
        iSize_t taskDim = mpi.m_lastTask - mpi.m_firstTask+1;
        dcmplx* temp = new dcmplx[taskDim];
        NonInteractingOflModel*  blochTable4 = blochTable;
        dcmplx* p_temp = temp;
        iSize_t counter = 0;
        //  Display a progress bar
        utilities::LoadBar progress;  
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            progress.Initialize(mpi.m_lastTask - mpi.m_firstTask+1);
        }
        for(kState_t k4=0; k4<kMax; ++k4, ++blochTable4)
        {   
            NonInteractingOflModel* blochTable3 = blochTable4;
            //  Make temporary tables for the spin-up/spin down Bloch eigenvectors
            blochTable4->GetBlochCoefficients(blochDup, 0, _SPIN_UP_, _NO_CONJUGATE_);
            blochTable4->GetBlochCoefficients(blochDdown, 0, _SPIN_DOWN_, _NO_CONJUGATE_);
            for(kState_t k3=k4; k3<kMax; ++k3, ++blochTable3)
            {
                NonInteractingOflModel*blochTable2 = blochTable;
                //  Make temporary tables for the spin-up/spin down Bloch eigenvectors
                blochTable3->GetBlochCoefficients(blochCup, 0, _SPIN_UP_, _NO_CONJUGATE_);
                blochTable3->GetBlochCoefficients(blochCdown, 0, _SPIN_DOWN_, _NO_CONJUGATE_);
                for(kState_t k2=0; k2<kMax; ++k2, ++blochTable2, ++counter)
                {
                    if(mpi.m_firstTask <= (int)counter && mpi.m_lastTask >= (int)counter)
                    {
                        blochTable2->GetBlochCoefficients(blochBup, 0, _SPIN_UP_, _CONJUGATE_);
                        blochTable2->GetBlochCoefficients(blochBdown, 0, _SPIN_DOWN_, _CONJUGATE_);
                        //  For every k2,k3,k4 there is a unique k1 in the BZ
                        //  such that k1+k2-k3-k4 = a reciprocal lattice vector
                        //  as G = gx G1 + gy G2
                        int gx;
                        int gy;
                        kState_t k1 = this->FindK1(k2, k3, k4, gx, gy);
                        blochTable[k1].GetBlochCoefficients(blochAup, 0, _SPIN_UP_, _CONJUGATE_);
                        blochTable[k1].GetBlochCoefficients(blochAdown, 0, _SPIN_DOWN_, _CONJUGATE_);
                        //  Evaluate the convolutions of the Bloch coefficients
                        const int offset1 = std::max(-gx,0);
                        const int offset2 = std::max(-gy,0);
                        const int offset3 = std::max(gx,0);
                        const int offset4 = std::max(gy,0);
                        utilities::ArrayConvolution2D<dcmplx, dcmplx, sizeof(dcmplx)>(blochBdown, blochDup, blochTableXDim,
                                                                                      blochTableYDim, convolution1, 0, 0, 0, 0);
                        utilities::ArrayConvolution2D<dcmplx, dcmplx, sizeof(dcmplx)>(blochCdown, blochAup, blochTableXDim, 
                                                                                      blochTableYDim, convolution2, offset1, 
                                                                                      offset2, offset3, offset4);
                        utilities::ArrayConvolution2D<dcmplx, dcmplx, sizeof(dcmplx)>(blochBup, blochDdown, blochTableXDim, 
                                                                                      blochTableYDim, convolution3, 0, 0, 0, 0);
                        utilities::ArrayConvolution2D<dcmplx, dcmplx, sizeof(dcmplx)>(blochCup, blochAdown, blochTableXDim, 
                                                                                      blochTableYDim, convolution4, offset1,
                                                                                      offset2, offset3, offset4);
                        //  Combine the convolutions (use Kahan accumulation to improve summation accuracy)
                        utilities::KahanAccumulation<dcmplx> Vkkkk;
                        for(iSize_t i=0; i<convolutionDim-offset1*convolutionArrayOffset-offset2 - offset3*convolutionArrayOffset-offset4; ++i)
                        {
                            Vkkkk += convolution1[i+offset3*convolutionArrayOffset+offset4]*convolution2[i+offset1*convolutionArrayOffset+offset2];
                            Vkkkk += convolution3[i+offset3*convolutionArrayOffset+offset4]*convolution4[i+offset1*convolutionArrayOffset+offset2];
                        }
                        //  Add the current term to the overall list of terms
                        *p_temp = Vkkkk.m_sum;
                        ++p_temp;
                        if(0 == mpi.m_id)	// FOR THE MASTER NODE
                        {
                            progress.Display(counter+1);
                        }
                    }
                }
            }
        }
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.AdditionalInfo()<<std::endl;
        }        
        //  clear up memory
        delete[] blochAup;
        delete[] blochAdown;
        delete[] blochBup;
        delete[] blochBdown;
        delete[] blochCup;
        delete[] blochCdown;
        delete[] blochDup;
        delete[] blochDdown;
        delete[] convolution1;
        delete[] convolution2;
        delete[] convolution3;
        delete[] convolution4;
	    MPI_Barrier(mpi.m_comm);
        MPI_Status status;
        mpi.Gather<dcmplx>(temp, taskDim, listOfQuarticTerms->data(), dimension, 0, mpi.m_comm, status);
        delete[] temp;
        return; 
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief This function converts the table format for storing Vkkkk, Ekk and 
    //! the momentum conserving tables from a regular array format to a multi
    //! key hash table format, if not done so already
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::ConvertTableFormat(
        const utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class 
    {
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t- CONVERTING LOOK-UP ARRAYS TO HASH TABLES"<<std::endl;
            }
            m_quadraticHashTables.Clear();
            m_quarticHashTables.Clear();
            m_quadraticHashTables.Initialize(m_params.m_dimX*m_params.m_dimY);
            m_quarticHashTables.Initialize(m_params.m_dimX*m_params.m_dimY);
            m_quadraticHashTables.SetFromArray(&m_quadraticTables, m_params.m_dimX*m_params.m_dimY);
            m_quarticHashTables.SetFromArray(&m_quarticTables, m_params.m_dimX*m_params.m_dimY);
            m_quarticTables.Clear();
            m_quadraticTables.Clear();
            m_params.m_tableFormat = _HASH_;
            m_quadraticHashTables.MpiSynchronize(0, mpi);
            m_quarticHashTables.MpiSynchronize(0, mpi);
        }
        return;
    }

    //!
    //! Store all tabulated coefficient data in a file
    //!
    void InteractingOflModel::TermsToFile(
        const io::fileFormat_t format,  //!<    Format of file
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<"\n\t- STORING QUARTIC AND QUADRATIC COEFFICIENT DATA IN FILES."<<std::endl;
        }
        std::stringstream fileNameQuadratic, filenameQuartic;
        fileNameQuadratic.str("");
        filenameQuartic.str("");
        fileNameQuadratic << m_params.m_outPath << "/quadratic_coefficient_table_kx_" << m_params.m_dimX << "_ky_";
        fileNameQuadratic << m_params.m_dimY << ".dat";
        filenameQuartic << m_params.m_outPath << "/quartic_coefficient_table_kx_" << m_params.m_dimX << "_ky_";
        filenameQuartic << m_params.m_dimY << ".dat";
        if(_ARRAY_ == m_params.m_tableFormat)
        {
            m_quadraticTables.ToFile(fileNameQuadratic.str(), format, mpi);
            m_quarticTables.ToFile(filenameQuartic.str(), format, mpi);
        }
        else
        {
            m_quadraticHashTables.ToFile(fileNameQuadratic.str(), format, mpi);
            m_quarticHashTables.ToFile(filenameQuartic.str(), format, mpi);
        }
        MPI_Barrier(mpi.m_comm);
	    return;
    }

    //!
    //! Retrieve all tabulated coefficient data from a file
    //! 
    void InteractingOflModel::TermsFromFile(
        const io::fileFormat_t format,  //!<    Format of file
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.SecondaryOutput()<<"\n\t- RETRIEVING QUARTIC AND QUADRATIC COEFFICIENT DATA FROM FILES...";
        }
        std::stringstream fileNameQuadratic, filenameQuartic;
        fileNameQuadratic.str("");
        filenameQuartic.str("");
        fileNameQuadratic << m_params.m_outPath << "/quadratic_coefficient_table_kx_" << m_params.m_dimX << "_ky_";
        fileNameQuadratic << m_params.m_dimY << ".dat";
        filenameQuartic << m_params.m_outPath << "/quartic_coefficient_table_kx_" << m_params.m_dimX << "_ky_";
        filenameQuartic << m_params.m_dimY << ".dat";
        if(_ARRAY_ == m_params.m_setTableFormat)
        {
            m_quarticTables.FromFile(filenameQuartic.str(), format, mpi);
            m_quadraticTables.FromFile(fileNameQuadratic.str(), format, mpi);
        }
        else
        {
            m_quarticHashTables.FromFile(filenameQuartic.str(), format, mpi);
            m_quadraticHashTables.FromFile(fileNameQuadratic.str(), format, mpi);
        }
        MPI_Barrier(mpi.m_comm);
        m_params.m_termTablesBuilt = true;
        m_params.m_tableFormat = m_params.m_setTableFormat;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<" - DONE"<<std::endl;
        }
	    return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Implement a change of basis to a cylindrical geometry using the basis
    //! of hybrid localized Wannier states as defined in Eq 32 of PRB 86,085129 (2012)
    //!
    //! This function updates the quartic (Vkkkk) and quadratic e_kk terms to a form 
    //! appropriate for the new basis (note that the quadratic term is no longer 
    //! diagonal in the new basis).
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
    //!
    //! The we substitute this change of basis into the existing Hamiltonian
    //! and calculate the new coefficients
    //!
    //! Since we no longer conserve kx momentum, the coefficient look-up tables
    //! are now more conveniently stored as a multi-key hash multi-maps
    //!
    //! Tables of Bloch coefficients used to calculate the factors phi, lambda
    //! and A_x should be passed as an argument to this function
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::ChangeToWannierBasis(
        NonInteractingOflModel* blochTable,  //!<    Lookup table of Bloch wave functions
        const double matrixElementTol,          //!<    Tolerance for setting matrix elements to precisely zero
        utilities::MpiWrapper& mpi)             //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.MainOutput()<<"\n\tRECALCULATING COEFFICIENTS FOR WANNIER BASIS CASE"<<std::endl;
        }
        //  Clear any existing hash tables, in preparation for population with new
        //  sets of coefficients
        m_quadraticHashTables.Clear();
        m_quarticHashTables.Clear();
        m_quadraticHashTables.Initialize(m_params.m_dimX*m_params.m_dimY);
        m_quarticHashTables.Initialize(m_params.m_dimX*m_params.m_dimY);
        //  Generate tables of angular momentum conserving values for the Wannier basis
        this->GenerateWannier2KTable(m_quadraticHashTables.GetKTable(), mpi);
        this->GenerateWannier4KTable(m_quarticHashTables.GetKTable(), mpi);
        m_params.m_tableFormat = _HASH_;

        //////      Recalculate quadratic terms     ////////////////////////////////////
        //  The new set of quadratic terms will be of the general form
        //  E_{x1,x2,ky}. It will be diagonal in ky still, but non-diagonal in x.
        //  sum_{kx,ky} E_{kx,ky} c^+_{kx,ky} c_{kx,ky}
        //  = sum_{x1,x2,ky} [sum_{kx} E_{kx,ky} |1/coefficient(kx,ky)|^2 1/dimX
        //  e^{2*I*PI*kx/dimX [x1-x2]}] c^+_{x1,ky} c_{x2,ky}
        //  The || of the coefficients is just 1 (since these coefficients are
        //  just phase factors. So we have:
        //
        //  E_{x1,x2,ky} = e^{2*I*PI*offsetX/dimX [x1-x2]} sum_{kx} E_{kx,ky} 1/dimX e^{2*I*PI*kx/dimX [x1-x2]}
        //
        //  The sums can be performed by a set of N_y discrete Fourier transforms
        //  The Fourier coefficients depend only on the difference x2-x1
        //  so we need only calculate the transform for all possible values
        //  of this difference (by Hermiticity we only need to consider x2>=x1).
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t- RECALCULATING QUADRATIC TERMS\n"<<std::endl;
            }
            utilities::LoadBar progress;
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                progress.Initialize(m_params.m_dimX*m_params.m_dimY);
            }
            //  Declare some variables to store input/outputs of Fourier transforms
            dcmplx* fourierCoefficients = new dcmplx[m_params.m_dimX*m_params.m_dimY];
            dcmplx* outputCoefficients  = new dcmplx[m_params.m_dimX*m_params.m_dimY];
            for(kState_t ky=0; ky<m_params.m_dimY; ++ky)
            {
                for(kState_t kx=0; kx<m_params.m_dimX; ++kx)
                {
                    fourierCoefficients[ky*m_params.m_dimX+kx] = m_quadraticTables.GetVTable()->at(kx*m_params.m_dimY+ky)/((double)m_params.m_dimX);
                    if(0 == mpi.m_id)	// FOR THE MASTER NODE
                    {
                        progress.Display(ky*m_params.m_dimX+kx+1);
                    }
                }
            }
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<std::endl;
            }
            utilities::DiscreteFourierTransform1D(m_params.m_dimX, m_params.m_dimY, fourierCoefficients, 
                                                  outputCoefficients, FFTW_BACKWARD);
            //  The output from the discrete transform is a set of coefficients like E_{ky,x2-x1}
            //  Next, we use that output to populate the look-up hash tables E_{x1,x2,ky}
            for(kState_t x1=0; x1<m_params.m_dimX; ++x1)
            {
                for(kState_t x2=0; x2<m_params.m_dimX; ++x2)
                {
                    for(kState_t ky=0; ky<m_params.m_dimY; ++ky)
                    {
                        m_quadraticHashTables.GetVTable()->Insert(utilities::Key(x1*m_params.m_dimY+ky, x2*m_params.m_dimY+ky)) = exp(2.0*I*PI*fabs(x2-x1)*m_params.m_offsetX/(double)m_params.m_dimX)*outputCoefficients[ky*m_params.m_dimX+abs(x2-x1)];
                    }
                }
            } 
            delete[] fourierCoefficients;
            delete[] outputCoefficients;
        }
        //////      Recalculate quartic terms     //////////////////////////////////////
        //  The new set of quadratic terms will be of the general form
        //
        //  V_{x1,ky1,x2,ky2,x3,ky3,x4,ky4}.
        //
        //  The list will contain a larger number of terms in the new basis
        //  because there is no conservation of the x values (so we need to 
        //  do the full sum over x1,x2,x3,x4). One of the ky values remains
        //  fixed by the other three, however, due to momentum conservation
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.SecondaryOutput()<<"\n\t- RECALCULATING QUARTIC TERMS\n"<<std::endl;
            }
            //////      Calculate and store the change of basis matrix
            NonInteractingOflModelGrid modelGrid;
            modelGrid.GenerateWannierCoefficients(m_params.m_dimX, m_params.m_dimY, m_params.m_offsetX, blochTable, 
                                                  &m_basisChangeMatrix[0], mpi);
            iSize_t dimBasisChange = m_params.m_dimX*m_params.m_dimY;
            //  Reserve a memory block with which to retrieve the momentum conserving list of k-values
            kState_t kRetrieveBuffer[m_quarticHashTables.GetMaxKCount()];
            //  The change of basis calculation is parallelized so that the 
            //  outermost sums (over x4 and k4y) are distributed over the available 
            //  nodes
            //  The hash tables are combined at the end of the calculation
            mpi.DivideTasks(mpi.m_id, m_params.m_dimX*m_params.m_dimY, mpi.m_nbrProcs,
                            &mpi.m_firstTask, &mpi.m_lastTask, false);
            iSize_t counter = 0;
            utilities::LoadBar progress;
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                progress.Initialize(mpi.m_lastTask+1);
            }
            for(kState_t x4=0; x4<m_params.m_dimX; ++x4)
            {
                for(kState_t k4y=0; k4y<m_params.m_dimY; ++k4y, ++counter)
                {
                    kState_t k4 = x4*m_params.m_dimY+k4y;
                    if(mpi.m_firstTask <= (int)counter && mpi.m_lastTask >= (int)counter)
                    {
                        if(0 == mpi.m_id)	// FOR THE MASTER NODE
                        {
                            progress.Display(counter+1);
                        }
                        for(kState_t x3=0; x3<m_params.m_dimX; ++x3)
                        {
                            for(kState_t k3y=0; k3y<m_params.m_dimY; ++k3y)
                            {
                                kState_t k3 = x3*m_params.m_dimY+k3y;
                                for(kState_t x2=0; x2<m_params.m_dimX; ++x2)
                                {
                                    for(kState_t k2y=0; k2y<m_params.m_dimY; ++k2y)
                                    {
                                        kState_t k2 = x2*m_params.m_dimY+k2y;
                                        //  Get the corresponding list of k1 values
                                        //  that satisfy y momentum conservation
                                        iSize_t nbrK1;
                                        m_quarticHashTables.GetK1(kRetrieveBuffer,nbrK1,k2,k3,k4);
                                        for(iSize_t j=0; j<nbrK1; ++j)
                                        {
                                            kState_t k1 = kRetrieveBuffer[j];
                                            kState_t k1y = utilities::Modulo<int>(k1, m_params.m_dimY);
                                            kState_t x1  = k1/m_params.m_dimY;      
                                            dcmplx summation = 0.0;     
                                            for(kState_t k4x=0; k4x<m_params.m_dimX; ++k4x)
                                            {
                                                for(kState_t k3x=0; k3x<m_params.m_dimX; ++k3x)
                                                {
                                                    for(kState_t k2x=0; k2x<m_params.m_dimX; ++k2x)
                                                    {
                                                        kState_t q1;
                                                        kState_t q2 = k2x*m_params.m_dimY+k2y;
                                                        kState_t q3 = k3x*m_params.m_dimY+k3y;
                                                        kState_t q4 = k4x*m_params.m_dimY+k4y;
                                                        //  Set the value of the q1 state using
                                                        //  kx,ky momentum conservation
                                                        //  (note that ky momentum conservation 
                                                        //  is already satisfied by k1y)
                                                        iSize_t nbrQ1;
                                                        m_quarticTables.GetK1(&q1, nbrQ1, q2, q3, q4);
                                                        //  Note that this value is unique
                                                        unsigned int uIndex1 = (x1*m_params.m_dimY+k1y)*dimBasisChange+q1;
                                                        unsigned int uIndex2 = (x2*m_params.m_dimY+k2y)*dimBasisChange+q2;
                                                        unsigned int uIndex3 = (x3*m_params.m_dimY+k3y)*dimBasisChange+q3;
                                                        unsigned int uIndex4 = (x4*m_params.m_dimY+k4y)*dimBasisChange+q4;
                                                        summation += m_quarticTables.GetVkkkk(q1, q2, q3, q4)*std::conj(m_basisChangeMatrix[uIndex1]*m_basisChangeMatrix[uIndex2])*m_basisChangeMatrix[uIndex3]*m_basisChangeMatrix[uIndex4];
                                                    }
                                                }
                                            }
                                            if(abs(summation)>matrixElementTol)
                                            {
                                                m_quarticHashTables.GetVTable()->Insert(utilities::Key(k1, k2, k3, k4)) = summation;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.AdditionalInfo()<<std::endl;
            }
        }
        //  Combine Hash tables populated on each node
        m_quarticHashTables.GetVTable()->MpiGather(0, mpi);
        m_quadraticHashTables.GetVTable()->MpiGather(0, mpi);
        m_quarticHashTables.GetVTable()->MpiSynchronize(0, mpi);
        m_quadraticHashTables.GetVTable()->MpiSynchronize(0, mpi);
        //  Clear existing non-hash tables
        m_quarticTables.Clear();
        m_quadraticTables.Clear();
        MPI_Barrier(mpi.m_comm);
    }
    
    //!
    //! Calculate the sigma^z map and store in m_magnetization
    //!
    void InteractingOflModel::CalcualteSigmaZMap(
        NonInteractingOflModel* blochTable,  //!<    Lookup table of Bloch wave functions
        const utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<"\n\t- GENERATING SIGMA^Z MAP";
            fflush(stdout);
        }
        NonInteractingOflModelGrid modelGrid;
        m_magnetization.resize(m_params.m_dimX*m_params.m_dimY);
        modelGrid.CalcualteSigmaZMap(m_params.m_dimX, m_params.m_dimY, blochTable, &m_magnetization[0]);
        m_params.m_magnetizationCalculated = true;
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<" - DONE"<<std::endl;
        }
        return;
    }

    //!
    //! Generate a Fock basis 
    //! 
    void InteractingOflModel::BuildFockBasis(
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    { 
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {   
            if(m_params.m_blockDiagonalize)
            {
                if(m_params.m_useWannierBasis)
                {
                    utilities::cout.MainOutput()<<"\n\t============ BUILDING FOCK BASIS (ky_tot = "<<m_params.m_kyTot<<") ============ "<<std::endl;
                }
                else
                {
                    utilities::cout.MainOutput()<<"\n\t============ BUILDING FOCK BASIS (kx_tot = "<<m_params.m_kxTot<<", ky_tot = "<<m_params.m_kyTot<<") ============ "<<std::endl;
                }
            }
            else
            {
                utilities::cout.MainOutput()<<"\n\t============ BUILDING FULL FOCK BASIS ============ "<<std::endl;
            }
        }
        //  Build the Fock basis first
        //  The function "GenerateFockSpace" requires a single argument:
        //  a function that takes a Fock basis state as an argument and returns 
        //  a bool saying whether that state is in the specified sector or not. 
        bool isZeroDimensional = false;
        if(m_params.m_blockDiagonalize)
        {
            //  To make the function we're going to use std::bind to fix parameters 
            //  in the function called TestLinearMomentumSector defined in
            //  linear_momentum_constraints_2D.h
            if(m_params.m_useWannierBasis)
            {
                //  Divide into ky sectors only
                FermionFockBasis basis;
                isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, m_params.m_dimX*m_params.m_dimY, 
                                                            std::bind(TestLinearMomentumSector1, std::placeholders::_1, 
                                                                      m_params.m_kyTot, m_params.m_dimY), mpi);
                m_hamiltonian.SetFockBasis(basis);
            }
            else
            {
                //  Divide into both kx and ky sectors
                FermionFockBasis basis;
                isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, m_params.m_dimX*m_params.m_dimY,
                                                            std::bind(TestLinearMomentumSector2, std::placeholders::_1, 
                                                                      m_params.m_kxTot, m_params.m_kyTot, m_params.m_dimX, 
                                                                      m_params.m_dimY), mpi);
                m_hamiltonian.SetFockBasis(basis);
            }
        }
        else
        {
            //  If we're not block diagonalizing then we pass a dummy function
            //  to m_hamiltonian (the function is declared as a "lambda function")
            //  This function will always return true.
            FermionFockBasis basis;
            isZeroDimensional = basis.GenerateFockSpace(m_params.m_nbrParticles, m_params.m_dimX*m_params.m_dimY, [](const fock_t state){return true;}, mpi);
            m_hamiltonian.SetFockBasis(basis);
        }
        if(isZeroDimensional)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            { 
                utilities::cout.SecondaryOutput()<<"\n\tERROR: FOCK BASIS IS ZERO DIMENSIONAL ON ONE OR MORE NODES "<<std::endl;
            }
            m_params.m_fockBasisBuilt = false;
            return;
        }
        m_params.m_fockBasisBuilt = true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a matrix representation of our Hamiltonian in the 2D momentum
    //! space basis.
    //!
    //! This implements the Hamiltonian given e.g. in eq (5) of PRL 109,265301.
    //!
    //! When called multiple times, this function will not construct a new Hamiltonian
    //! unless the previous Hamiltonian has been cleared with the ClearHamiltonian
    //! function.
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::BuildHamiltonian(
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        if(m_params.m_termTablesBuilt && m_params.m_fockBasisBuilt && !m_params.m_hamiltonianBuilt)
        {
            //////  GENERATE THE HAMILTONIAN        //////
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {   
                utilities::cout.MainOutput()<<"\n\t============ BUILDING HAMILTONIAN ============ "<<std::endl;
            }
            m_hamiltonian.Initialize(m_params.m_nbrParticles, m_params.m_dimX*m_params.m_dimY, utilities::_SPARSE_MAPPED_, mpi);
            MPI_Barrier(mpi.m_comm);
            //////      Build the Hamiltonian by including different types of term
            if(_ARRAY_ == m_params.m_tableFormat)
            {
                m_hamiltonian.Add_CdC_Terms(&m_quadraticTables, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                m_hamiltonian.Add_CdCdCC_Terms(&m_quarticTables, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
            }
            else
            {
                m_hamiltonian.Add_CdC_Terms(&m_quadraticHashTables, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
                m_hamiltonian.Add_CdCdCC_Terms(&m_quarticHashTables, 0, m_hamiltonian.m_data.m_fockSpaceDim, mpi);
            }
            m_hamiltonian.PrintMemoryAllocation(mpi);
            for(int i=0; i<mpi.m_nbrProcs; ++i)
            {
                if(mpi.m_id == i)
                {
                    m_hamiltonian.PrintHamiltonian(i, mpi);
                }
                MPI_Barrier(mpi.m_comm);
            }
            MPI_Barrier(mpi.m_comm);
            m_params.m_hamiltonianBuilt = true;
        }
        else
        {
            m_params.m_hamiltonianBuilt = false;
        }
    }

    //!
    //! Diagonalize our Hamiltonian using a selected method
    //! 
    void InteractingOflModel::Diagonalize(
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {
        if(m_params.m_hamiltonianBuilt)
        {
            //////  CALL SELECTED DIAGONALIZATION ROUTINE
            //  Switch to LAPACK for very small matrices
            const unsigned int smallMatrixLimit = 1000;
            if(m_hamiltonian.m_data.m_fockSpaceDim<smallMatrixLimit)
            {
                if(0 == mpi.m_id)
                {
                    utilities::cout.SecondaryOutput()<<"\n\tMATRIX DIMENSION <"<<smallMatrixLimit<<" - SETTING DIAGONALIZATION METHOD TO FULL\n"<<std::endl;
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
                    if(m_params.m_useWannierBasis)
                    {
                        inFileName<<"_sector_"<<m_params.m_kyTot;
                        outFileName<<"_sector_"<<m_params.m_kyTot;
                    }
                    else
                    {
                        inFileName<<"_sector_"<<m_params.m_kxTot<<"_"<<m_params.m_kyTot;
                        outFileName<<"_sector_"<<m_params.m_kxTot<<"_"<<m_params.m_kyTot;
                    }
                }
                
                inFileName<<"_proc_"<<mpi.m_id<<".bin";
                outFileName<<"_proc_"<<mpi.m_id<<".bin";
                m_hamiltonian.UpdateArpackFileNames(inFileName.str(),outFileName.str(),mpi);
                if(_ARRAY_ == m_params.m_tableFormat)
                {
                    MatrixVectorFunction<dcmplx, QuadraticTermTablesBase<dcmplx>, QuarticTermTables<dcmplx>> mv(m_hamiltonian, mpi);
                    m_params.m_hamiltonianDiagonalized = m_hamiltonian.LanczosDiagonalize(
                    &m_quadraticTables, &m_quarticTables, mv, m_params.m_nbrEigenvaluesToFind, mpi);
                }
                else
                {
                    MatrixVectorFunction<dcmplx, QuadraticTermHashTablesBase<dcmplx>, QuarticTermHashTablesBase<dcmplx>> mv(m_hamiltonian, mpi);
                    m_params.m_hamiltonianDiagonalized = m_hamiltonian.LanczosDiagonalize(
                    &m_quadraticHashTables, &m_quarticHashTables, mv, m_params.m_nbrEigenvaluesToFind, mpi);
                }
            }
            if(!m_params.m_hamiltonianDiagonalized)
            {
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::cerr << "\n\tWARNING - NOT ABLE TO RUN DIAGONALIZATION ROUTINE." << std::endl;
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Call a function to plot the Hamiltonian
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::PlotHamiltonian(
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
    {       
        return m_hamiltonian.PlotHamiltonian(m_params.m_outFileName, m_params.m_outPath, mpi);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Output all stored eigenvalue and eigenvector data to a text file
    //!
    //! This function MUST be called on all nodes if run in parallel 
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::EigensystemToFile(
        const bool writeEigenvalues,    //!<    Option to write eigenvalues to file
        const bool writeEigenvectors,   //!<    Option to write eigenvectors to file
        const io::fileFormat_t format,  //!<    Format of file
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
            ss << m_params.MakeBaseFileName(io::_OUT_) << "_eigensystem.dat";
            m_hamiltonian.EigensystemToFile(ss.str(), writeEigenvalues, writeEigenvectors, format, mpi);
            mpi.ExitFlagTest();
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Obtain eigenvalue and eigenvector data from existing files
    //!
    //! This function MUST be called on all nodes if run in parallel 
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::EigensystemFromFile(
        const bool readEigenvalues,         //!<    Option to read eigenvalues from file
        const bool readEigenvectors,        //!<    Option to read eigenvalues from file
        const io::fileFormat_t format,      //!<    Format of file
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
    {
        if((readEigenvalues || readEigenvectors) && !m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {   
                utilities::cout.SecondaryOutput()<<"\n\t============ READING EIGENSYSTEM DATA FROM FILE ============ "<<std::endl;
            }
            m_hamiltonian.Initialize(m_params.m_nbrParticles, m_params.m_dimX*m_params.m_dimY, utilities::_SPARSE_MAPPED_, mpi);
            //  If the matrix is Block-diagonalized then generate an appropriate set
            //  of allowed Fock basis states
            if(m_params.m_blockDiagonalize)
            {
                //  For an explanation of this function call see the comments
                //  in the BuildHamiltonian() function
                FermionFockBasis basis;
                basis.GenerateFockSpace(m_params.m_nbrParticles, m_params.m_dimX*m_params.m_dimY, std::bind(TestLinearMomentumSector2, 
                                        std::placeholders::_1, m_params.m_kxTot, m_params.m_kyTot, 
                                        m_params.m_dimX, m_params.m_dimY), mpi);
                m_hamiltonian.SetFockBasis(basis);
            }
            std::stringstream ss;
            ss.str("");
            ss << m_params.MakeBaseFileName(io::_IN_) << "_eigensystem.dat";
            m_hamiltonian.EigensystemFromFile(ss.str(), readEigenvalues, readEigenvectors, format, mpi);
            mpi.ExitFlagTest();
            //  Counts as if the Hamiltonian were diagonalized
            m_params.m_hamiltonianDiagonalized = true;
        }
        return;
    }

    //!
    //! Store the full Hamiltonian in (a) file(s) - one file for each node
    //! that the storage is distributed over. CRS format is used 
    //!
    void InteractingOflModel::HamiltonianToFile(
        const io::fileFormat_t format,  //!<    Format of file 
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        std::stringstream fileName;
        fileName.str("");
        fileName << m_params.MakeBaseFileName(io::_OUT_) << "_proc_" << mpi.m_id;
        fileName << "_hamiltonian.dat";
        m_hamiltonian.HamiltonianToFile(fileName.str(), format, mpi);
        return;
    }

    //!
    //! Retrieve the full Hamiltonian from file(s) - one file for each node
    //! that the storage is distributed over. CCS format is used 
    //!
    void InteractingOflModel::HamiltonianFromFile(
        const io::fileFormat_t format,  //!<    Format of file
        utilities::MpiWrapper& mpi)     //!<    Instance of the mpi wrapper class
    {
        m_hamiltonian.Initialize(m_params.m_nbrParticles,m_params.m_dimX*m_params.m_dimY,utilities::_SPARSE_CRS_,mpi);
        std::stringstream fileName;  
        fileName.str("");    
        fileName << m_params.MakeBaseFileName(io::_IN_) << "_proc_" << mpi.m_id;
        fileName << "_hamiltonian.dat";
        m_hamiltonian.HamiltonianFromFile(fileName.str(), format, mpi);
        m_params.m_hamiltonianBuilt = true;
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to clear any currently stored FermionHamiltonian object and
    //! reset the m_params.m_hamiltonianDiagonalized and m_params.m_hamiltonianBuilt flags
    ////////////////////////////////////////////////////////////////////////////////
    void InteractingOflModel::ClearHamiltonian()
    {
        if(m_params.m_hamiltonianBuilt)
        {
            m_hamiltonian.Clear();
            m_params.m_hamiltonianBuilt = false;
            m_params.m_hamiltonianDiagonalized = false;
        }
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief A function to determine the total probability for occupying each quantum
    //! state corresponding to a specified eigenvector of the interacting Hamiltonian.
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::CalculateOccupations(
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.MainOutput()<<"\n\tGENERATING OCCUPATION DATA\n"<<std::endl;
            }
            double* probabilityList = new double[m_params.m_dimX*m_params.m_dimY];
            for(iSize_t eigenvectorIndex=0;eigenvectorIndex<m_params.m_nbrEigenvaluesToFind;++eigenvectorIndex)
            {
                //////      GATHER DATA      /////////////////////////////////////////
                //  Call a function in the FermionHamiltonian class to bin the occupation
                //  probabilities
                Observables observables;
                observables.DensityFunction(m_hamiltonian, eigenvectorIndex, probabilityList, mpi);
                MPI_Barrier(mpi.m_comm);
                //////      WRITE PROBABILITY DATA TO FILE
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::ofstream f_out;
                    std::stringstream fileName;
                    fileName.str("");
                    fileName<<m_params.MakeBaseFileName(io::_OUT_);
                    fileName<<"_eig_"<<eigenvectorIndex<<"_occupations.dat";
                    f_out.open(fileName.str().c_str(), std::ios::out);
                    if(!f_out.is_open())
                    {
                        std::cerr<<"\n\tERROR in CalculateOccupations: Cannot open file "<<fileName.str()<<std::endl;
                        mpi.m_exitFlag = true;
                    }
                    else
                    {
                        //  Map out the unit cell
                        utilities::LatticeVector2D<double> kVector(m_params.m_offsetX, m_params.m_offsetY);
                        utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1, (double)1/m_params.m_dimX);
                        utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2, (double)1/m_params.m_dimY);
                        f_out.precision(15);
                        f_out<<m_params.m_dimX<<"\n";
                        f_out<<m_params.m_dimY<<"\n";
                        double* p_prob = probabilityList;
                        for(iSize_t x=0; x<m_params.m_dimX; ++x, kVector-=reciprocalLattice::G2,kVector += increment1)
                        {
                            for(iSize_t y=0; y<m_params.m_dimY; ++y, kVector += increment2, ++p_prob)
                            {
                                f_out<<std::fixed<< kVector.GetKx()<<"\t"<<kVector.GetKy()<<"\t"<<*(p_prob)<<"\n";
                            }
                        }
                        f_out.close();
                        utilities::cout.SecondaryOutput()<<"\n\tOccupations data output to file "<<fileName.str()<<std::endl;
                    }
                }
                mpi.ExitFlagTest();
            }
            delete[] probabilityList; 
            return true;
        }
        else
        {
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                std::cerr<<"\n\t- CalculateOccupations ERROR: eigenvalues/eigenvectors not allocated/no longer allocted!"<<std::endl;
            }
            return false;
        }
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the value of the mean squared magnetization
    //!
    //! sum_{k,k'} n_k n_k' sigma_k sigma_k' 
    //!
    //! Note that since the mean magnetization is zero, this value will also 
    //! give the magnetic susceptibility of the model.
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::CalculateSusceptibility(
        utilities::MpiWrapper& mpi)       //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized && m_params.m_magnetizationCalculated)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.MainOutput()<<"\n\tGENERATING SUSCEPTIBILITY DATA\n"<<std::endl;
            }
            iSize_t nbrOrbitals = m_hamiltonian.m_data.m_nbrOrbitals;
            double* probabilityList = new double[nbrOrbitals*nbrOrbitals];
            for(iSize_t eigenvectorIndex=0; eigenvectorIndex<m_params.m_nbrEigenvaluesToFind; ++eigenvectorIndex)
            {
                //////      GATHER DATA      /////////////////////////////////////////
                //  Call a function in the FermionHamiltonian class to evaluate
                //  an array of of n_k n_k'
                Observables observables;
                observables.DensityDensityFunction(m_hamiltonian, eigenvectorIndex, probabilityList, mpi);
                MPI_Barrier(mpi.m_comm);
                //  Calculate a weighted sum of n_k n_k' sigma_k sigma_k'
                double totalSquaredMagnetization = 0.0;
                for(iSize_t i=0; i<nbrOrbitals; ++i)
                {
                    for(iSize_t j=0; j<nbrOrbitals; ++j)
                    {
                        totalSquaredMagnetization += probabilityList[i*nbrOrbitals+j]*m_magnetization[i]*m_magnetization[j];
                    }
                }
                //////      Write susceptibility data to file
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::ofstream f_out;
                    std::stringstream fileName;
                    fileName.str("");
                    fileName<<m_params.MakeBaseFileName(io::_OUT_);
                    fileName<<"_eig_"<<eigenvectorIndex<<"_susceptibility.dat";
                    f_out.open(fileName.str().c_str(), std::ios::out);
                    if(!f_out.is_open())
                    {
                        std::cerr<<"\n\tERROR in CalculateSusceptibility: Cannot open file "<<fileName.str()<<std::endl;
                        mpi.m_exitFlag = true;
                    }
                    else
                    {
                        f_out<<std::fixed<<totalSquaredMagnetization<<"\n";
                        f_out.close();
                        utilities::cout.SecondaryOutput()<<"\n\tSusceptibility data output to file "<<fileName.str()<<std::endl;
                    }
                }
                mpi.ExitFlagTest();
            }
            delete[] probabilityList;
            return true;
        }
        else
        {
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                if(!m_params.m_hamiltonianDiagonalized)
                {
                    std::cerr<<"\n\t- CalculateSusceptibility ERROR: eigenvalues/eigenvectors not allocated/no longer allocted!"<<std::endl;
                }
                else
                {
                    std::cerr<<"\n\t- CalculateSusceptibility ERROR: Magnetization not calculated!"<<std::endl;
                }
            }
            return false;
        }
        return false;
    }
                
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate a list of the nbrStates most probable Fock states. The
    //! ||^2 amplitude associated with the states is output along with their 
    //! Fock space index. 
    //!
    //! The output file contains a numerical representation of the state
    //! (which should be interpreted as a binary occupation number)
    //! followed by the real and imaginary parts of the amplitude). 
    //!
    //! e.g. For a 4x4 system with n=3 we might have:
    //!
    //! ## Most probable state in eigenvector 0
    //! 17472           -0.486148      -0.161189
    //!
    //! 17472 here is to be interpreted as the occupation pattern:
    //! 0100,0100,0100,0000
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::ListMostProbable(
        const iSize_t nbrStates,    //!<    Number of the most probable states to output
        utilities::MpiWrapper& mpi)         //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)    //  For the master node
            {
                utilities::cout.MainOutput()<<"\n\t============ LISTING MOST PROBABLE FOCK STATES ============"<<std::endl;
            }
            //  Allocate space to store the most probable amplitudes associated with 
            //  each eigenvector
            fock_t* mostProbableStates = 0;
            dcmplx* mostProbableAmplitudes  = 0;
            if(0 == mpi.m_id)    //  For the master node
            {
                mostProbableStates      = new fock_t[nbrStates*m_params.m_nbrEigenvaluesToFind];
                mostProbableAmplitudes  = new dcmplx[nbrStates*m_params.m_nbrEigenvaluesToFind];
            }
            //  Extract the most probable states from the eigenvector data
            Observables observables; 
            observables.GetMostProbable(m_hamiltonian, nbrStates, mostProbableStates, mostProbableAmplitudes, mpi);
            if(0 == mpi.m_id)    //  For the master node
            {
                //  Prepare to output a list of squared overlap values into a file
                std::stringstream outFileName;
                outFileName.str("");
                outFileName<<m_params.MakeBaseFileName(io::_OUT_)<<"_most_probable.dat";
                std::ofstream f_probable;
                f_probable.open(outFileName.str().c_str(), std::ios::out);
                if(!f_probable.is_open())
                {
                    std::cerr<<"\n\tERROR in ListMostProbable: Cannot open file "<<outFileName.str()<<std::endl;
                    mpi.m_exitFlag = true;
                }
                else
                {
                    utilities::cout.SecondaryOutput()<<"\n\t - Writing most probable state data to file "<<outFileName.str()<<std::endl;
                    //  Write the data to a text file in a format consistent with 
                    //  reading into the CalculateOverlaps function
                    for(iSize_t e=0; e<m_params.m_nbrEigenvaluesToFind; ++e)
                    {
                        //  Generate a text description of the state
                        for(iSize_t j=0; j<nbrStates; ++j)
                        {
                            f_probable<<"## ";
                            switch(j)
                            {
                                case 0:
                                    f_probable<<"Most";
                                    break;
                                case 1:
                                    f_probable<<"2nd most";
                                    break;
                                case 2:
                                    f_probable<<"3rd most";
                                    break;
                                default:
                                    f_probable<<j+1<<"th most";
                                    break;
                            }
                            f_probable<<" probable state in eigenvector "<<e<<"\n";
                            //  Output the Fock state and its amplitude
                            f_probable<<std::setw(15)<<std::left<<mostProbableStates[j+e*nbrStates];
                            f_probable<<"\t"<<std::setw(15)<<std::left<<std::real(mostProbableAmplitudes[j+e*nbrStates]);
                            f_probable<<std::setw(15)<<std::left<<std::imag(mostProbableAmplitudes[j+e*nbrStates])<<"\n\n";
                        }
                    }
                    f_probable.close();
                }
            } 
            if(0 != mostProbableStates)     delete[] mostProbableStates;
            if(0 != mostProbableAmplitudes) delete[] mostProbableAmplitudes;
            mpi.ExitFlagTest();
            return true;
        }
        return false;
    }
                 
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the density-density function and store it in a file
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::CalculateDensityDensityFunction(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.MainOutput()<<"\n\tGENERATING DENSITY-DENSITY FUNCTION DATA\n"<<std::endl;
            }
            iSize_t nbrOrbitals = m_hamiltonian.m_data.m_nbrOrbitals;
            double* probabilityList = new double[nbrOrbitals*nbrOrbitals];
            for(iSize_t eigenvectorIndex=0; eigenvectorIndex<m_params.m_nbrEigenvaluesToFind; ++eigenvectorIndex)
            {
                //////      GATHER DATA      /////////////////////////////////////////
                //  Call a function in the FermionHamiltonian class to evaluate
                //  an array of of n_k n_k'
                Observables observables;
                observables.DensityDensityFunction(m_hamiltonian, eigenvectorIndex, probabilityList, mpi);
                MPI_Barrier(mpi.m_comm);
                //////      Write density-density function data to file
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::ofstream f_out;
                    std::stringstream fileName;
                    fileName.str("");
                    fileName<<m_params.MakeBaseFileName(io::_OUT_);
                    fileName<<"_eig_"<<eigenvectorIndex<<"_density_density.dat";
                    f_out.open(fileName.str().c_str(), std::ios::out);
                    if(!f_out.is_open())
                    {
                        std::cerr<<"\n\tERROR in CalculateDensityDensityFunction: Cannot open file "<<fileName.str()<<std::endl;
                        mpi.m_exitFlag = true;
                    }
                    else
                    {
                        //  Map out the unit cell for the k and k' indices
                        utilities::LatticeVector2D<double> kVector(m_params.m_offsetX, m_params.m_offsetY);
                        utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1, (double)1/m_params.m_dimX);
                        utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2, (double)1/m_params.m_dimY);
                        f_out.precision(15);
                        f_out<<m_params.m_dimX<<"\n";
                        f_out<<m_params.m_dimY<<"\n";
                        f_out<<std::setw(20)<<std::left<<"kx1 values"<<"\t"<<std::setw(20)<<std::left<<"ky1 values"<<std::setw(20)<<std::left<<"kx2 values"<<"\t"<<std::setw(20)<<std::left<<"ky2 values"<<"\t"<<std::setw(20)<<std::left<<"density-density"<<"\n";
                        double* p_prob = probabilityList;
                        for(iSize_t x=0; x<m_params.m_dimX; ++x, kVector-=reciprocalLattice::G2, kVector += increment1)
                        {
                            for(iSize_t y=0; y<m_params.m_dimY; ++y, kVector += increment2)
                            {
                                utilities::LatticeVector2D<double> qVector(m_params.m_offsetX, m_params.m_offsetY);
                                for(iSize_t x1=0; x1<m_params.m_dimX; ++x1, qVector-=reciprocalLattice::G2, qVector += increment1)
                                {
                                    for(iSize_t y1=0; y1<m_params.m_dimY; ++y1, qVector += increment2, ++p_prob)
                                    {
                                        f_out<<std::fixed<<kVector.GetKx()<<"\t"<<kVector.GetKy()<<"\t"<<qVector.GetKx()<<"\t"<<qVector.GetKy()<<"\t"<<*(p_prob)<<"\n";
                                    }
                                }
                            }
                        }
                        f_out.close();
                        utilities::cout.SecondaryOutput()<<"\n\tDensity-Density function data output to file "<<fileName.str()<<std::endl;
                    }
                }
                mpi.ExitFlagTest();
            }
            delete[] probabilityList;  
            return true;
        }
        else
        {
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                std::cerr<<"\n\t- CalculateDensityDensityFunction ERROR: eigenvalues/eigenvectors not allocated/no longer allocted!"<<std::endl;
            }
            return false;
        }
        return false;
    }
               
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the inverse participation ratio and store it in a file
    //!
    //! Inverse participation ratio is defined as 1/ sum |psi|^4.
    //! This function gives a measure of the number of orbitals over
    //! which the many body wave function may be considered to be
    //! localized
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::CalculateParticipationRatio(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.MainOutput()<<"\n\tGENERATING PARTICIPATION RATIO DATA\n"<<std::endl;
            }
            
            for(iSize_t eigenvectorIndex=0; eigenvectorIndex<m_params.m_nbrEigenvaluesToFind; ++eigenvectorIndex)
            {
                Observables observables;
                double ratio = observables.GetParticipationRatio(m_hamiltonian, eigenvectorIndex, mpi);
                MPI_Barrier(mpi.m_comm);
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    //  Prepare a file to be written to
                    std::ofstream f_out;
                    std::stringstream fileName;
                    fileName.str("");
                    fileName<<m_params.MakeBaseFileName(io::_OUT_);
                    fileName<<"_eig_"<<eigenvectorIndex<<"_participation_ratio.dat";
                    f_out.open(fileName.str().c_str(), std::ios::out);
                    if(!f_out.is_open())
                    {
                        std::cerr<<"\n\tERROR in CalculateParticipationRatio: Cannot open file "<<fileName.str()<<std::endl;
                        mpi.m_exitFlag = true;
                    }
                    else
                    {                        
                        f_out.precision(15);
                        f_out<<ratio<<"\n";
                        f_out.close();
                        utilities::cout.SecondaryOutput()<<"\n\tParticipation ratio data output to file "<<fileName.str()<<std::endl;
                    }
                }
                mpi.ExitFlagTest();
            }
            return true;
        }
        return false;
    }
                 
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Calculate the modified k-space correlation function
    //!
    //! The function is given by:
    //!
    //! f(G) = sum_{k1,k2} < c^+_k2-G c_k2 c^+_k1+G c_k1 >
    //!
    //! This function attempts to measure the existence of special density
    //! wave configurations, which might show up as peaks at particular values
    //! in reciprocal space, G. G covers 4x4 unit cells in this implementation
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::CalculateTranslationalDensityDensityFunction(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.MainOutput()<<"\n\tGENERATING TRANSLATIONAL DENSITY-DENSITY FUNCTION DATA\n"<<std::endl;
            }
            iSize_t nbrOrbitals = m_hamiltonian.m_data.m_nbrOrbitals;
            const iSize_t nbrUnitCells = 1;
            std::vector<double> functionList(nbrUnitCells*nbrUnitCells*nbrOrbitals);
            for(iSize_t eigenvectorIndex=0;eigenvectorIndex<m_params.m_nbrEigenvaluesToFind;++eigenvectorIndex)
            {
                //////      GATHER DATA      /////////////////////////////////////////
                //  Call a function in the FermionHamiltonian class to evaluate
                //  sum_{k1,k2} < c^+_k2-G c_k2 c^+_k1+G c_k1 >
                Observables observables;
                observables.TranslationalDensityDensityFunction(m_hamiltonian, eigenvectorIndex,
                    std::bind(StateAddition, std::placeholders::_1, std::placeholders::_2, m_params.m_dimX, 
                              m_params.m_dimY, nbrUnitCells*m_params.m_dimX, nbrUnitCells*m_params.m_dimY),
                    std::bind(StateSubtraction, std::placeholders::_1, std::placeholders::_2, m_params.m_dimX, 
                              m_params.m_dimY, nbrUnitCells*m_params.m_dimX, nbrUnitCells*m_params.m_dimY),
                    functionList, mpi);
                MPI_Barrier(mpi.m_comm);
                //////      Write function to a file
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::ofstream f_out;
                    std::stringstream fileName;
                    fileName.str("");
                    fileName<<m_params.MakeBaseFileName(io::_OUT_);                   
                    fileName<<"_eig_"<<eigenvectorIndex<<"_translational_density_density.dat";
                    f_out.open(fileName.str().c_str(), std::ios::out);
                    if(!f_out.is_open())
                    {
                        std::cerr<<"\n\tERROR in CalculateTranslationalDensityDensityFunction: Cannot open file "<<fileName.str()<<std::endl;
                        mpi.m_exitFlag = true;
                    }
                    else
                    {
                        //  Map out the unit cell
                        utilities::LatticeVector2D<double> kVector(m_params.m_offsetX, m_params.m_offsetY);
                        utilities::LatticeVector2D<double> increment1(reciprocalLattice::G1, (double)1/m_params.m_dimX);
                        utilities::LatticeVector2D<double> increment2(reciprocalLattice::G2, (double)1/m_params.m_dimY);
                        f_out.precision(15);
                        f_out<<nbrUnitCells*m_params.m_dimX<<"\n";
                        f_out<<nbrUnitCells*m_params.m_dimY<<"\n";
                        auto it_func = functionList.begin();
                        for(iSize_t x=0; x<nbrUnitCells*m_params.m_dimX; ++x, kVector-=reciprocalLattice::G2, kVector += increment1)
                        {
                            for(iSize_t y=0; y<m_params.m_dimY*nbrUnitCells; ++y, kVector += increment2, ++it_func)
                            {
                                f_out<<std::fixed<< kVector.GetKx()<<"\t"<<kVector.GetKy()<<"\t"<<*(it_func)<<"\n";
                            }
                        }
                        f_out.close();
                        utilities::cout.SecondaryOutput()<<"\n\tTranslational density-density data output to file "<<fileName.str()<<std::endl;
                    }  
                }
                mpi.ExitFlagTest();
            }
            return true;
        }
        else
        {
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                std::cerr<<"\n\t- CalculateTranslationalDensityDensityFunction ERROR: eigenvalues/eigenvectors not allocated/no longer allocated!"<<std::endl;
            }
            return false;
        }
        return false;
    }
                
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief Generate a list of the total rotational density-density operator 
    //! associated with  a pair of orbitals for a given eigenvector
    //!
    //! Given by:
    //!
    //! R(m) =  1/6N sum_{n,k} exp(i*m*n*pi/3) <c^+_{R_n k} c_{R_n k} c^+_k c_k>
    //!
    //! Where n,m in {0,1,2,3,4,5} and R is a map of orbitals to their 60 degree 
    //! rotated positions
    //!
    //! \return true if the calculation was successful, false otherwise
    ////////////////////////////////////////////////////////////////////////////////
    bool InteractingOflModel::CalculateRotationalDensityDensityFunction(
        utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
        const
    {
        if(m_params.m_hamiltonianDiagonalized)
        {
            if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                utilities::cout.MainOutput()<<"\n\tGENERATING ROTATIONAL DENSITY-DENSITY FUNCTION DATA\n"<<std::endl;
            }
            dcmplx functionList[6];
            for(iSize_t eigenvectorIndex=0; eigenvectorIndex<m_params.m_nbrEigenvaluesToFind; ++eigenvectorIndex)
            {
                //////      GATHER DATA      /////////////////////////////////////////
                //  Call a function in the FermionHamiltonian class to evaluate R(m)
                Observables observables;
                observables.RotationalDensityDensityFunction(m_hamiltonian, eigenvectorIndex,
                    std::bind(R60, std::placeholders::_1, std::placeholders::_2, m_params.m_dimX, m_params.m_dimY),
                    std::bind(InverseR60, std::placeholders::_1, std::placeholders::_2, m_params.m_dimX, m_params.m_dimY),
                    functionList, mpi);
                MPI_Barrier(mpi.m_comm);
                //////      Write function to a file
                if(0 == mpi.m_id)	// FOR THE MASTER NODE
                {
                    std::ofstream f_out;
                    std::stringstream fileName;
                    fileName.str("");
                    fileName<<m_params.MakeBaseFileName(io::_OUT_);
                    fileName<<"_eig_"<<eigenvectorIndex<<"_rotational_density_density.dat";
                    f_out.open(fileName.str().c_str(), std::ios::out);
                    if(!f_out.is_open())
                    {
                        std::cerr<<"\n\tERROR in CalculateRotationalDensityDensityFunction: Cannot open file "<<fileName.str()<<std::endl;
                        mpi.m_exitFlag = true;
                    }
                    else
                    {
                        for(int m=0; m<6; ++m)
                        {
                            f_out<<std::fixed<<m<<"\t"<<"\t"<<std::abs(functionList[m])<<"\n";
                        }
                        f_out.close();
                        utilities::cout.SecondaryOutput()<<"\n\tRotational density-density data output to file "<<fileName.str()<<std::endl;
                    }
                }
                mpi.ExitFlagTest();
            }
            return true;
        }
        else
        {
            if(0 == mpi.m_id) //  FOR THE MASTER NODE
            {
                std::cerr<<"\n\t- CalculateRotationalDensityDensityFunction ERROR: eigenvalues/eigenvectors not allocated/no longer allocated!"<<std::endl;
            }
            return false;
        }
        return false;
    }
}   //  End namespace diagonalization
