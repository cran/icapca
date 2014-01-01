/* Copyright 1995-2013 Roger P. Woods, M.D. */
/* Modified 7/4/13 */

/*
 * compute the minimum and maximum bounds for the argument in the
 * gamma function
 *
 * xmin and xmax are not the only bounds, but they are the only non-
 * trivial ones to calculate
 *
 * On exit
 *	*xmin=minimum legal value of x in gamma(x). Any smaller values
 *		might lead to underflow
 *	*xmax=maximum legal value of x in gamma(x)...overflow
 */

#include "AIR.h"
#include <float.h>

AIR_Error AIR_dgamlm(double *xmin, double *xmax)

{
	{
		double alnsml=log(DBL_MIN);
		*xmin=-alnsml;
		{
			unsigned int i;
			
			for(i=1;i<=10;i++){
			
				double xold=*xmin;
				double xln=log(*xmin);
				
				*xmin=*xmin-*xmin*((*xmin+.5)*xln-*xmin-.2258+alnsml)/(*xmin*xln+.5);
				if(fabs(*xmin-xold)<.005) break;
			}
			if(i>10){
#ifndef USING_R
				printf("%s: %d: ",__FILE__,__LINE__);
				printf("unable to find xmin\n");
#else
				REprintf("%s: %d: ",__FILE__,__LINE__);
				REprintf("unable to find xmin\n");
#endif
				return AIR_DGAMLM_CANT_MIN_ERROR;
			}
		}
		*xmin=-*xmin+.01;
	}
	{
		double alnbig=log(DBL_MAX);
		*xmax=alnbig;
		{
			unsigned int i;
			
			for(i=1;i<=10;i++){
			
				double xold=*xmax;
				double xln=log(*xmax);
				
				*xmax=*xmax-*xmax*((*xmax-.5)*xln-*xmax+.9189-alnbig)/(*xmax*xln-.5);
				if(fabs(*xmax-xold)<.005) break;
			}
			if(i>10){
#ifndef USING_R
				printf("%s: %d: ",__FILE__,__LINE__);
				printf("unable to find xmax\n");
#else
				REprintf("%s: %d: ",__FILE__,__LINE__);
				REprintf("unable to find xmax\n");
#endif
				return AIR_DGAMLM_CANT_MAX_ERROR;
			}
		}
		*xmax=*xmax-.01;
	}
	if(-*xmax+1.0>*xmin) *xmin=-*xmax+1.0;
	return 0;
}

	
