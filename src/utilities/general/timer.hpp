////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 21/03/2014
//!
//!  \file
//!		This file contains functions to record the time taken to perform 
//!     a task and print this information
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

#ifndef _SIMON_TIMER_HPP_INCLUDED_
#define _SIMON_TIMER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "cout_tools.hpp" //  Functions to manipulate std::cout output
#include <chrono>       //  A c++11 library for formatting time displays
#include <ctime>        //  For strftime

namespace utilities
{
    ////////////////////////////////////////////////////////////////////////////////
    //! \brief  A class to provide a simple timer to record the time between the 
    //! class being constructed and destroyed.
    //!
    ////////////////////////////////////////////////////////////////////////////////

    class Timer
    {
        private:
        
        std::chrono::high_resolution_clock::time_point m_wallTime;
        std::clock_t m_cpuTime;
    
        ////////////////////////////////////////////////////////////////////////////////
        //! \brief    TimeStamp prints the current YMDHMS date
        //!
        ////////////////////////////////////////////////////////////////////////////////

        inline void TimeStamp()
        {
	        static char time_buffer[40];
	        const struct tm *tm_ptr;
	        time_t now;

	        now = time(NULL);
	        tm_ptr = localtime ( &now );

	        strftime(time_buffer,40,"%d %B %Y %I:%M:%S %p",tm_ptr);

	        utilities::cout.AdditionalInfo() << time_buffer << "\n\n";

	        return;
        }
        
        public:

		//!
		//! Start timing function
		//!
		inline void Start()
		{
		    utilities::cout.AdditionalInfo()<<"\n\t TIMING COMMENCED ";
		    TimeStamp();
		    m_cpuTime  = std::clock();
            m_wallTime = std::chrono::high_resolution_clock::now();
		}
    
        //!
        //! Step timing function
        //!
        inline void Stop()
        {
            auto duration = std::chrono::high_resolution_clock::now() - m_wallTime;
            auto seconds  = std::chrono::duration_cast<std::chrono::seconds>(duration);
            auto millis   = std::chrono::duration_cast<std::chrono::milliseconds>(duration-seconds);

            utilities::cout.AdditionalInfo()<<"\n\t TIMING ENDED ";
		    TimeStamp();
		    utilities::cout.AdditionalInfo()<<"\t CPU TIME ELAPSED "<<(double)(std::clock() - m_cpuTime) / CLOCKS_PER_SEC<<" SECONDS.\n"<<std::endl;
            utilities::cout.AdditionalInfo()<<"\t WALL TIME ELAPSED "<<seconds.count()<<"."<<millis.count()<<" SECONDS.\n"<<std::endl;
           
		}
    };
}

//  End namespace utilities

#endif

