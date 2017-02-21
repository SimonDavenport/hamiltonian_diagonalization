////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport 
//!                                                                             
//!	 \file	
//!     This header file contains a static constant definition of sqrt(-1)
//!     as a complex double type
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef _I_CONST_DEFINED_
#define _I_CONST_DEFINED_
#include <complex>
static const std::complex<double> I=std::complex<double> (0, 1);
#endif
