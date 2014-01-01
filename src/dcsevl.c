/* Copyright 1995-2013 Roger P. Woods, M.D. */
/* Modified 7/4/13 */

/* double dcsevl()
 *
 * Evaluate the N-term Chebyshev series CS at X.
 * From the Fortran program in LINPACK library by the same name
 *
 * Input
 *	X	value at which the series is to be evaluated
 *	CS	array of N terms of a Chebyshev series. In evaluating
 *		 CS, only half the first cofficient is summed
 *	N	number of terms in array CS.
 *
 * Returns:
 *	The value of interest if no error
 *	errcode if error
 */

#include "AIR.h"
#include <float.h>

double AIR_dcsevl(const double x, const double *cs, const unsigned int n, AIR_Error *errcode)

{
	static double 	onepl;
	static AIR_Boolean first=TRUE;
	
	/*Load constants*/
	if(first){
		onepl=1.0+DBL_EPSILON;
		first=FALSE;
	}
	if(n==0){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("function called with fewer than one term\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("function called with fewer than one term\n");
#endif
		*errcode=AIR_DCSEVL_BAD_N_ERROR;
		return 0.0;
	}
	if(n>1000){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("function called with more than 1000 terms\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("function called with more than 1000 terms\n");
#endif
		*errcode=AIR_DCSEVL_BAD_N_ERROR;
		return 0.0;
	}
	if(fabs(x)>onepl){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("x outside the interval (-1,+1)\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("x outside the interval (-1,+1)\n");
#endif
		*errcode=AIR_DCSEVL_BAD_X_ERROR;
		return 0.0;
	}
	{
		double b1=0.0;
		double b0=0.0;
		double twox=2.0*x;
		double b2;
		
		{
			unsigned int i=n-1;
			
			do{
				b2=b1;
				b1=b0;
				b0=twox*b1-b2+cs[i];
			}while(i--!=0);
		}
		*errcode=0;
		return 0.5*(b0-b2);
	}
}

	
