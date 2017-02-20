////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 03/08/2014
//!
//!  \file
//!		This file contains a class that accepts a string as the contents of a
//!     script. The class contains a function to write the string to a temporary
//!     file name, tagged uniquely by the system time, a function to attempt to
//!     execute the script, and a function to remove the file name after execution
//!     is completed. 
//!
//!                    Copyright (C) 2014 Simon C Davenport
//!
//!		This program is free software: you can redistribute it and/or modify
//!		it under the terms of the GNU General Public License as published by
//!		the Free Software Foundation, either version 3 of the License,
//!		or (at your option) any later version.
//!
//!		This program is distributed in the hope that it will be useful, but
//!		WITHOUT ANY WARRANTY; without even the implied warranty of 
//!		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
//!		General Public License for more details.
//!
//!		You should have received a copy of the GNU General Public License
//!		along with this program. If not, see <http://www.gnu.org/licenses/>.
//! 
////////////////////////////////////////////////////////////////////////////////

#ifndef _RUN_SCRIPT_HPP_INCLUDED_
#define _RUN_SCRIPT_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <string>       //  For std::string
#include <fstream>      //  For file output
#include <sstream>      //  For std::stringstream
#include <ctime>        //  For obtaining the system time
#include "cout_tools.hpp" //  For manipulating cout

namespace utilities
{

//////////////////////////////////////////////////////////////////////////////////
//! \brief A class to convert a given string to a text file and execute it as a 
//! script. 
//!
//! The class contains a function to write the string to a temporary file name, 
//! tagged uniquely by the system time, a function to attempt to execute the 
//! script, and a function to remove the file name after execution is completed.
//!
//////////////////////////////////////////////////////////////////////////////////

class Script
{
    private:

    std::string m_scriptText;
    std::string m_fileNameStem;
    std::string m_tempFileName;

    void GenFileName()
    {
        //  Tag the script file name with the full system time
    
        std::stringstream tempFileName;
        
        tempFileName.str("");

        tempFileName.precision(10);
        
        time_t now;
        
        time(&now);
        
        tempFileName<<m_fileNameStem<<"_"<<now<<".py";
        
        m_tempFileName = tempFileName.str();
    };
    
    inline void RunSystemCommand(const std::string command) const
    {
        int sysReturn = system(command.c_str());
        
        utilities::cout.DebuggingInfo()<<"\n\tScript run info: system returned value "<<sysReturn<<std::endl;
    };
    
    inline void ScriptToFile() const
    {
        //  Put the script in an executable output file
    
        std::ofstream f_script;
 
        f_script.open(m_tempFileName.c_str(),std::ios::out);

        if(!f_script.is_open())
        {
            std::cerr<<"ERROR: file "<<m_tempFileName<<" could not be generated"<<std::endl;
            exit(EXIT_FAILURE);
        }

        f_script<<m_scriptText.c_str();
        
        f_script.close();
        
        //  Run a system command to change the executable bit on the output file
        
        std::stringstream command;
        
        command.str("");
        
        command<<"chmod +x "<<m_tempFileName;
        
        this->RunSystemCommand(command.str());
    };

    public:
    
    //  Default constructor
    
    Script()
    :
        m_scriptText(""),
        m_fileNameStem("./temp_script_file")
    {
        this->GenFileName();
    };
    
    inline void SetScript(const std::string newText)
    {
        m_scriptText = newText;
    };
    
    inline void SetFileNameStem(const std::string newStem)
    {
        std::stringstream stem;
        
        stem.str("");
        
        stem<<"./"<<newStem;
        
        m_fileNameStem = stem.str();
        
        this->GenFileName();
    };
    
    inline void Execute() const
    {
        //  put the current script text in an executable file
    
        this->ScriptToFile();
        
        //  call the system command to run the script

        this->RunSystemCommand(m_tempFileName);
        
        //  call the system command to remove the temporary file
        
        std::stringstream command;
        
        command.str("");
        
        command<<"rm "<<m_tempFileName;
        
        this->RunSystemCommand(command.str());
    };

};  //  End namespace utilities

};

#endif
