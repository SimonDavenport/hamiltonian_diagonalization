////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 19/06/2014
//!
//!  \file
//!		This file contains a wrapper for the FFTW library
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "fftw_wrapper.hpp"

//////////////////////////////////////////////////////////////////////////////// 
//! \brief A namespace to contain any functions and utilities that I have 
//! written for use with any c++ program.
//!
////////////////////////////////////////////////////////////////////////////////

namespace utilities
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//	

    //////////////////////////////////////////////////////////////////////////////////
	//!	\brief This function performs a series of 1D discrete Fourier transforms
	//!
	//! Y_{j,m} = sum_n=1^N exp^{sign*2*I*PI*n*j/N} X_{n,m}
	//!
	//////////////////////////////////////////////////////////////////////////////////

    void DiscreteFourierTransform1D(
        int N,          //!<    Dimension of the sum in the transform
        int M,          //!<    Number of independent transforms performed in one function call
        dcmplx* inputs, //!<    Array of Fourier coefficients for each transform
                        //!     Must be of dimension N*M
        dcmplx* outputs,//!<    Array to contain the output data. Must be of
                        //!     dimension N*M
        int sign)       //!<    Set the sign of the transform as either
                        //!     FFTW_BACKWARD (+1) or FFTW_FORWARD (-1)
    {
        //  declare the input/output arrays and the fftw_plan object
        
        fftw_complex *in, *out;
        fftw_plan plan = 0;
        
        //  Allocate memory for input and output arrays
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
        
        for(int m=0;m<M;m++)
        {
            //  Prepare the "plan" (use the "FFTW_PATIENT" method)

            plan = fftw_plan_dft_1d(N,in,out,sign,FFTW_PATIENT);
        
            //  Set the input data
            
            for(int n=0;n<N;n++)
            {
                in[n][0] = std::real(inputs[m*N+n]);
                in[n][1] = std::imag(inputs[m*N+n]);
            }
            
            //  Execute the plan
            
            fftw_execute(plan);

            //  Extract the output
            
            for(int n=0;n<N;n++)
            {
                outputs[m*N+n] = dcmplx(out[n][0],out[n][1]);
            }
        }
        
        //  Clear up memory allocations
        
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
        
        return;
    }

};  //  End namespace utilities 

