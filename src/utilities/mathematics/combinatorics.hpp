////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!  	Header file for implementations of binomial and factorial functions,
//!     and sets of permutations.	
//!
//!                    Copyright (C) Simon C Davenport
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

#ifndef _COMBINATORICS_HPP_INCLUDED_
#define _COMBINATORICS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <vector>
#include <algorithm>
#include <iostream>
#if _ENABLE_HIGH_PRECISION_
#include "../wrappers/high_precision_wrapper.hpp"
#endif

namespace utilities
{
    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function returns the factorial of the input
    //!	
    //! Tempalte parameter T specifies calcualtion/output type (e.g. using
    //! long int will cause overflow errors around 20!, so it can be set
    //! to use doubles instead)
    //!
    //!	\return	 Value of factorial
    //////////////////////////////////////////////////////////////////////////////////
    template<typename T>
    T Factorial(
	    const int input)    //!<	The number you want to take the factorial of
    {
	    if(input == 0)
	    {
		    return 1;
	    }
	    T output = 1;
	    for(int k=1; k<=input; ++k)
	    {
		    output *= k;
	    }
	    return output;
    }

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_HIGH_PRECISION_

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief This function returns the ratio of the factorials.
	//!	To save processing time, common factors are cancelled out.
	//!	
	//!	The output argument in this case is an arbitrary precision variable wrapper
	//!	class HpWrap<type,precision> with integer 'precision' P)
	//!	
	//!	e.g. DivFactorial<mpft_t,128>(input1,input2,*output)
	//!	
	//!	Template argument U corresponds to the arbitrary precision type.
	//!
	//!	Template argument P corresponds to the precision level. 
	////////////////////////////////////////////////////////////////////////////////
	template <typename U, int P>
	void DivFactorial(
		double input1, 			//!<	factorial of input1 goes in the numerator
		double input2,			//!<	factorial of input2 goes in the denominator
		HpWrap<U,P>* output)	//!<	pointer to the output arbitrary precision wrapper object
	{
		double cmp;
		if(input1==input2)
		{
			(*output).Set(1.0);
			return;
		}
		HpWrap<U,P> tmp;
		cmp=1;
		if(input1>input2)
		{
			(*output).Set(input1);
			while (input1-cmp>input2)
			{
				(*output).MulBy((input1-cmp),*output);
				++cmp;
			}
		}
		else if(input2>input1)
		{
			(*output).Set(input2);
			while (input2-cmp>input1)
			{
				(*output).MulBy((input2-cmp),(*output));
				++cmp;
			}
			tmp.Set(1.0);
			tmp/=(*output);
			(*output)=tmp;
		}
		return;
	}
#endif
////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	//!	\brief An overload of the DivFactorial template function to avoid the
	//!	arbitrary  precision wrapper class.	
	//!
	//!	This function returns the ratio of the factorials.
	//!	To save processing time, common factors are cancelled out.
	//!	
	//!	e.g. DivFactorial<double>(input1,input2,*output)
	//!
	//!	\return Output of type U.
	////////////////////////////////////////////////////////////////////////////////
	template <typename T>
	T DivFactorial(
		T input1, 			//!<	factorial of input1 goes in the numerator
		T input2)			//!<	factorial of input2 goes in the denominator
	{
		T cmp;
		T output=0;
		if(input1==input2){output=1.0; return output;}
		cmp=1;
        if(input1<0)
        {
            return 0;
        }
		if(input1>input2)
		{
			output=input1;
			while(input1-cmp>input2)
			{
			    output*=(input1-cmp);
			    ++cmp;
			}
		}
		else if(input2>input1)
		{
			output=input2;
			while(input2-cmp>input1)
			{
			    output*=(input2-cmp);
			    ++cmp;
			}
			output=1/output;
		}
		return output;
	}

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This function returns the combinatoric factor input1 Choose input2.
    //!	Inputs are double variables to avoid overflows.
    //!	
    //!	\return	 Value of binomial coefficient (output type T)
    //////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    T Binomial(
	    T input1, 	//!<	First index of binomial
	    T input2)	//!<	Second index of binomial
    {
	    if(input2>input1||input2<0){return 0;}
	    if(input2==0)
	    {
		    return 1;
	    }
	    T cmp = 1;
	    T output = 0;
	    if(input1==input2){return 1;}
	    if(input1>input2){output=input1;
		    while (input1-cmp>input2){output*=(input1-cmp);++cmp;}}
	    else if(input2>input1){output=input2;
		    while (input2-cmp>input1){output*=(input2-cmp);++cmp;}
		    output=1 / output;}
	    T in=input1-input2;
	    cmp=0;
	    while (in-cmp>1){output/=(in-cmp);++cmp;}
	    return (T)round(output);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //! \brief	This template function generates the unique permutations of a list of 
    //! objects of arbitrary type T
    //////////////////////////////////////////////////////////////////////////////////
    template <class T>
    void UniqueObjectPermutations(
        std::vector<std::vector<T> >& objectPermutations,    
                                    //!<    List of all object permutations
                                    //!     to be populated
	    std::vector<T> objects)     //!<	List of objects to be permuted
    {
        std::sort(objects.begin(), objects.end());
	    do
	    {
            objectPermutations.push_back(objects);
        } 
        while (std::next_permutation(objects.begin(), objects.end()));
    }
}   //  End namespace utilities
#endif
