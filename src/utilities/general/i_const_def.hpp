////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport 
//!                                                                             
//!                      \date Last Modified: 14/02/2014
//!                                                                             
//!	 \file	
//!     This header file contains a static constant definition of sqrt(-1)
//!     as a complex double type
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef _I_CONST_DEFINED_
#define _I_CONST_DEFINED_

#include <complex>

static const std::complex<double> I=std::complex<double> (0,1);		
    //!<	Define the square root of -1

#endif
