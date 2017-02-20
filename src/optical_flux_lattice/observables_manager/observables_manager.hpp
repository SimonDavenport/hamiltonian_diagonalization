////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport
//!                                                                             
//!                      \date Last Modified: 27/01/2015
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

#ifndef _OBSERVABLES_MANAGER_HPP_INCLUDED_
#define _OBSERVABLES_MANAGER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <string>                       //  For std::string
#include <vector>                       //  For std::vector
#include <boost/program_options.hpp>    //  For program option handling
#include <map>                          //  For list of observables
#include "../optical_flux_lattice_hamiltonian/optical_flux_lattice_hamiltonian.hpp"

#include "../../utilities/wrappers/mpi_wrapper.hpp"   
                                        //  Include mpi functionality
#include "../../utilities/general/cout_tools.hpp"     
                                        //  Include cout functions

//  Include associated command line options
#include "../program_options/observables_options.hpp"

namespace diagonalization
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

class ObservablesManager
{
    private:
        
        std::map<std::string,bool> m_observablesList;
                                            //!<  A map specifying whether each observable
                                            //!   name should be measured or not   
        
        //////      Parameters associated with specific observables
        
        iSize_t m_nbrMostProbable;          //!<  The number of most probable values to store 
                                            //!   with the most-probable-list option
        iSize_t m_nbrA;                     //!<  Number of particles in the PES cut
        iSize_t m_orbitalCut;               //!<  Number of orbitals in the OES cut
        bool m_calcualtedObservables;       //!<  Flag set to true once
                                            //!   CalculateAllObservables is called
      
    public:

        ObservablesManager();
        
        ObservablesManager(ObservablesManager& other);

        void MpiSynchronize(const int nodeId,const utilities::MpiWrapper& mpi);

        void Print(const utilities::MpiWrapper& mpi);

        void AddObservables(boost::program_options::variables_map* optionList,utilities::MpiWrapper& mpi);
        
        int GetNbrSetObservables() const;
        
        void CalculateAllObservables(OpticalFluxLatticeHamiltonian* hamiltonian,utilities::MpiWrapper& mpi);
        
        void UpdateSqlFlags(boost::program_options::variables_map* optionList,const utilities::MpiWrapper& mpi);
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

}   //  End namespace diagonalization

#endif
