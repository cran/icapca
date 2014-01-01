/* Copyright 1995-2013 Roger P. Woods, M.D. */
/* Modified 7/4/13 */

/*
 * Initializes the orthogonal series, represented by the array OS, so
 * that INITDS is the number of terms needed to insure the error is
 * no larger than ETA. Ordinarily, ETA will be chosen to be one-
 * tenth machine precision.
 *
 * Modified from the FORTRANÊprogram of the same name in the
 * SLATEK library by W. Fullerton
 *
 * Input:
 *	OS 	double precision array of NOS coefficients in an
 *	 	orthogonal series.
 *
 *	NOS	number of coefficients in OS.
 *	ETA	single precision scalar containing requested
 *		accuracy of series.
 *
 * Returns:	The number of terms needed
 *		0, if error (errors are all WARNINGS)
 *
 */

#include "AIR.h"

unsigned int AIR_initds(const double *os, const unsigned int nos, const double eta)

{
	if(nos==0){
#ifndef USING_R
		printf("WARNING initds: called with less than one coefficient\n");
#else
		Rprintf("WARNING initds: called with less than one coefficient\n");
#endif
		return 0;
	}
	{
		unsigned int i;
		
		{
			double err=0.0;
			
			for(i=nos;(i--)!=0;){
			
				err+=fabs(os[i]);
				if(err>eta) break;
			}
		}
		if(i==nos-1){
#ifndef USING_R
			printf("WARNING initds: Chebyshev series too short for specified accuracy/n");
#else
			REprintf("WARNING initds: Chebyshev series too short for specified accuracy/n");
#endif
			return 0;
		}
		return i+1;
	}
}
	
