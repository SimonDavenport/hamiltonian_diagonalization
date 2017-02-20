////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 06/05/2015
//!
//!  \file
//!		This header file implements an array convolution of two arrays. 
//!     The function implemented here is identical to Mathematica's 
//!     "ListCorrelate" function, with maximal left and right overhang.
//!
//!     The Kahan Summation algorithm is implemented in order to reduce
//!     numerical errors.
//!
//!                    Copyright (C) 2015 Simon C Davenport
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

#ifndef _ARRAY_CONVOLUTION_HPP_INCLUDED
#define _ARRAY_CONVOLUTION_HPP_INCLUDED

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "kahan_arithmetic.hpp"         //  For Kahan accumulator class
#include "../general/template_tools.hpp"//  For template type checking
                                        //  funciton is_same
#include <iostream>                     //  For std::cerr

#if _ENABLE_HIGH_PRECISION_
#include "../wrappers/high_precision_wrapper.hpp"
                                        //  Enable use of arbitrary precision
                                        //  arithmetic
#endif

namespace utilities
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief A function to calculate and return the size of the convolved array,
    //! given the lengths of the input arrays and treating them as 1D arrays.
    //!
    //////////////////////////////////////////////////////////////////////////////// 

    inline int GetArrayConvolution1DSize(
        const int leftDim,          //!<    Dimension of left array
        const int rightDim)         //!<    Dimension of right array
    {
        return leftDim + rightDim - 1;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//    

    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief A template function to emulate Mathematica's "ListCorrelate" function
    //! for arbitrary c++ data types stored in regular 1D arrays.
    //!
    //! IMPORTANT: the output array MUST be allocated with dimension given by the 
    //! GetArrayConvolution1DSize function.
    //!
    //! This function is equivalent to the Mathematica function:
    //!
    //! ListCorrelate[rightArray,leftArray,{rightDim-startOffset,endOffset-rightDim},0]
    //!
    //! (assuming that rightDim<lefDim here)
    //!
    //////////////////////////////////////////////////////////////////////////////// 
    
    template <typename T>
    void ArrayConvolution1D(
        T* leftArray,           //!<    Left 1D array to be convolved
        const int leftDim,      //!<    Dimension of left array
        T* rightArray,          //!<    Right 1D array to be convolved
        const int rightDim,     //!<    Dimension of right array
        T* output,              //!<    Pointer to output array
        const int startOffset,  //!<    Skip the first startOffset terms
        const int endOffset)    //!<    Skip the final startOffset terms
    {
        T* p_output = output;
        
        const int convDim = GetArrayConvolution1DSize(leftDim,rightDim);
    
        for(int i=0;i<startOffset;++i,++p_output)
        {
            *p_output = 0;
        }
    
        for(int i=startOffset;i<convDim-endOffset;++i,++p_output)
        {
            //  Initialize to zero
            //*p_output = 0;
            
            KahanAccumulation<T> accum;
            
            //std::cout<<"\toutput = "<<*p_output;
            
            //  Dimension of the overlap
            
            int kernalDim = std::min(i+1,convDim-i);
            
            //  If the arrays are not of equal length, then restrict the 
            //  kernal of the convolution to take this into account
            
            if(kernalDim>leftDim)  kernalDim = leftDim;
            if(kernalDim>rightDim) kernalDim = rightDim;

            //std::cout<<"\ti="<<i<<" kernalDim = "<<kernalDim<<std::endl;
        
            if(rightDim>i)
            { 
                for(int j=0;j<kernalDim;j++)
                {
                    //std::cout<<"\tj="<<j<<std::endl;
                    //std::cout<<"\t+="<<leftArray[j]<<" * "<<rightArray[j+rightDim-kernalDim]<<std::endl;
                
                    //*p_output += leftArray[j]*rightArray[j+rightDim-kernalDim];
                    //std::cout<<"\toutput = "<<*p_output;
                    
                    accum += leftArray[j]*rightArray[j+rightDim-kernalDim];
                }
            }
            else
            {
                for(int j=0;j<kernalDim;j++)
                {
                    //std::cout<<"\tj="<<j+i-rightDim+1<<std::endl;
                    //std::cout<<"\t+="<<leftArray[j+i-rightDim+1]<<" * "<<rightArray[j]<<std::endl;
                
                    //*p_output += leftArray[j+i-leftDim+1]*rightArray[j];
                    //std::cout<<"\toutput = "<<*p_output;
                    
                    accum += leftArray[j+i-rightDim+1]*rightArray[j];
                }
            }

            #if 0
            
            if(abs(*p_output - accum.m_sum) > 0.000000000000001*abs(accum.m_sum))
            {
                std::cout<<" A DIFFERENCE?"<<std::endl;
                std::cout.precision(16);
                
                std::cout<<" ACTUAL "<<*p_output<<" VS ACCUM "<<accum.m_sum<<std::endl;
                
                getchar();
            }
            
            #endif
            
            *p_output = accum.m_sum;
        }
        
        for(int i=0;i<endOffset;++i,++p_output)
        {
            *p_output = 0;
        }
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//    
  
    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief A function to calculate and return the size of the convolved array,
    //! given the lengths of the input arrays and treating them as 2D arrays.
    //!
    //////////////////////////////////////////////////////////////////////////////// 

    inline int GetArrayConvolution2DSize(
        const int dimX,             //!<    Array x dimension
        const int dimY)             //!<    Array y dimension
    {
        return (2*dimY - 1 )*(2*dimX - 1);
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
   
    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief A template function to emulate Mathematica's "ListCorrelate" function
    //! for arbitrary c++ data types stored in regular arrays. This version of 
    //! the function treats the arrays as 2 dimensional.
    //!
    //! IMPORTANT: the output array MUST be allocated with dimension given by the 
    //! GetArrayConvolution2DSize function.
    //!
    //! This function is equivalent to the Mathematica function:
    //!
    //! ListCorrelate[rightArray,leftArray,{{dimX-topOffset,dimY-leftOffset},
    //! {bottomOffset-dimX,rightOffset-dimY}},0]
    //!
    //! where rightArray and leftArray are now indexed as 2D arrays. i.e.
    //! leftArray = {{a,b,c},{d,e,f}} etc.
    //!
    //! Standard arithmetic in convolution array calcualtion is improved through 
    //! the use of Kahan accumulation techniques. Alternatively, high precision 
    //! variables can be set using the template parameters.
    //!
    //! Template parameters: T - the type used in the convolved arrays and output
    //!                      H - a high precision mpc_t or mpfr_t type can be set 
    //!                      to be used in internal arithmetic. Otherwise T and H
    //!                      should be the same type
    //!                      P - an integer specifying the precision level
    //!
    //////////////////////////////////////////////////////////////////////////////// 
    
    template <typename T,typename H,int P>
    void ArrayConvolution2D(
        T* rightArray,          //!<    Right 1D array to be convolved
        T* leftArray,           //!<    Left 1D array to be convolved
        const int dimX,         //!<    Array x dimension
        const int dimY,         //!<    Array y dimension
                                //!     (for simplicity assuming both arrays have the same dimensions)
        T* output,              //!<    Pointer to output array
        const int topOffset,    //!<    Skip the first topOffset terms in each row
        const int leftOffset,   //!<    Skip the first leftOffset terms in each column
        const int bottomOffset, //!<    Skip the last bottomOffset terms in each row
        const int rightOffset)  //!<    Skip the final rightOffset terms in each column
    {
        //std::cout<<"\ttopOffset: "<<topOffset<<" leftOffset: "<<leftOffset<<" bottomOffset: "<<bottomOffset<<" rightOffset: "<<rightOffset<<std::endl;
    
        //std::cout<<"RIGHT ARRAY: "<<std::endl;
        //for(int i=0;i<dimX*dimY;i++)
        //{
        //    std::cout<<"rightArray["<<i<<"] = "<<rightArray[i]<<std::endl;
        //}
        
        //std::cout<<"LEFT ARRAY: "<<std::endl;
        //for(int i=0;i<dimX*dimY;i++)
        //{
        //    std::cout<<"leftArray["<<i<<"] = "<<leftArray[i]<<std::endl;
        //}
    
        T* p_output = output;
        
        const int convDimX = (2*dimX - 1);
        const int convDimY = (2*dimY - 1);
        
        //  Otherwise, evaluate the convolution
    
        for(int i=0;i<convDimX;i++)
        {
            //  X-dimension of the overlap        
            const int kernalDimX = std::min(i+1,convDimX-i);
        
            for(int j=0;j<convDimY;j++,p_output++)
            {   
                if(is_same<T,H>::value)
                {
                    //  Use standard precision arithmetic,
                    //  but make use of Kahan accumulation
                    
                    KahanAccumulation<T> accum;
                    
                    //  If the value is outside of the range specified by the offset parameters
                    //  then set the result to be zero
                    
                    if(i>= topOffset && i< convDimX - bottomOffset && j>=leftOffset && j< convDimY-rightOffset)
                    {
                        //  Y-dimension of the overlap      
                        const int kernalDimY = std::min(j+1,convDimY-j);
                        
                        //std::cout<<"\ti="<<i<<"\tj="<<j<<std::endl;
                        //std::cout<<"\tkernalDimX="<<kernalDimX<<"\tkernalDimY="<<kernalDimY<<std::endl;
                        
                        if(dimX>i)
                        { 
                            if(dimY>j)
                            {            
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        accum += leftArray[x*dimY+y]*rightArray[(x+dimX-kernalDimX)*dimY+(y+dimY-kernalDimY)];
                                    }
                                }
                            }
                            else
                            {
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        accum += leftArray[x*dimY+(y+dimY-kernalDimY)]*rightArray[(x+dimX-kernalDimX)*dimY+y];
                                    }
                                }
                            }
                        }
                        else
                        {
                            if(dimY>j)
                            {            
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        accum += leftArray[(x+dimX-kernalDimX)*dimY+y]*rightArray[x*dimY+(y+dimY-kernalDimY)];
                                    }
                                }
                            }
                            else
                            {
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        accum += leftArray[(x+dimX-kernalDimX)*dimY+(y+dimY-kernalDimY)]*rightArray[x*dimY+y];
                                    }
                                }
                            }
                        }
                    }
                    
                    *p_output = accum.m_sum;
                }
                #if _ENABLE_HIGH_PRECISION_
                else if((is_same<T,dcmplx>::value && is_same<H,mpfCmplx>::value) || (is_same<T,double>::value && is_same<H,mpf_t>::value))
                {
                    //  Determine that high precision complex variables are used with dcmplx in/out
                    //  or that high precision real variables are used with double in/out
                    
                    HpWrap<H,P> hpAccum;
                
                    hpAccum.Set(0.0);
                    
                    HpWrap<H,P> temp;
                    HpWrap<H,P> temp1;
                    
                    //  If the value is outside of the range specified by the offset parameters
                    //  then set the result to be zero
                    
                    if(i>= topOffset && i< convDimX - bottomOffset && j>=leftOffset && j< convDimY-rightOffset )
                    {
                        //  Y-dimension of the overlap      
                        const int kernalDimY = std::min(j+1,convDimY-j);
                        
                        //std::cout<<"\ti="<<i<<"\tj="<<j<<std::endl;
                        //std::cout<<"\tkernalDimX="<<kernalDimX<<"\tkernalDimY="<<kernalDimY<<std::endl;
                        
                        if(dimX>i)
                        { 
                            if(dimY>j)
                            {            
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        temp.Set(leftArray[x*dimY+y]);
                                        temp1.Set(rightArray[(x+dimX-kernalDimX)*dimY+(y+dimY-kernalDimY)]);
                                        temp *= temp1;
                                        hpAccum += temp;
                                    }
                                }
                            }
                            else
                            {
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        temp.Set(leftArray[x*dimY+(y+dimY-kernalDimY)]);
                                        temp1.Set(rightArray[(x+dimX-kernalDimX)*dimY+y]);
                                        temp *= temp1;
                                        hpAccum += temp;
                                    }
                                }
                            }
                        }
                        else
                        {
                            if(dimY>j)
                            {            
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        temp.Set(leftArray[(x+dimX-kernalDimX)*dimY+y]);
                                        temp1.Set(rightArray[x*dimY+(y+dimY-kernalDimY)]);
                                        temp *= temp1;
                                        hpAccum += temp;
                                    }
                                }
                            }
                            else
                            {
                                for(int x=0;x<kernalDimX;x++)
                                {
                                    for(int y=0;y<kernalDimY;y++)
                                    {
                                        temp.Set(leftArray[(x+dimX-kernalDimX)*dimY+(y+dimY-kernalDimY)]);
                                        temp1.Set(rightArray[x*dimY+y]);
                                        temp *= temp1;
                                        hpAccum += temp;
                                    }
                                }
                            }
                        }
                    }
                    
                    *p_output = hpAccum.Get();
                }
                #endif
                else
                {
                    std::cerr<<"\n\tERROR WITH ArrayConvolution2D - Template parameters T and H inconsistent"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                
                //std::cout<<"\toutput = "<<*p_output<<std::endl;getchar();
                #if 0
                if(abs(hpAccum.retrieve() - accum.m_sum) > 0.000000000000001*abs(accum.m_sum))
                {
                    std::cout<<" A DIFFERENCE?"<<std::endl;
                    
                    std::cout.precision(16);
                    
                    std::cout<<" ACTUAL "<<hpAccum.retrieve()<<" VS ACCUM "<<accum.m_sum<<std::endl;
                    
                    getchar();
                }
                
                #endif
            }    
        }
        
        return;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
}   //  End namespace utilities

#endif
