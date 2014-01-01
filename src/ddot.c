/* Copyright 1995-2002 Roger P. Woods, M.D. */
/* Modified 7/15/02 */

/*
 * double ddot()
 *
 * This routine will compute dot product of two vectors
 */

#include "AIR.h"

double AIR_ddot(const unsigned int n, const double *x, const double *y)

{
	/* Returns 0 if n==0 */
	double temp=0.0;
	{
		unsigned int i;

		for(i=0;i<n;i++){
			temp+=*x++**y++;
		}
	}
	return temp;
}
