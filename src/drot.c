/* Copyright 1996-2001 Roger P. Woods, M.D. */
/* Modified 5/12/01 */

/*
 * void drot()
 */


#include "AIR.h"

void AIR_drot(const unsigned int n, double *dx, double *dy, const double dc, const double ds)

{
	if(n==0 || (ds==0.0 && dc==1.0)) return;
	{
		unsigned int i;

		for(i=0;i<n;i++){
			double w=*dx;
			double z=*dy;
			*dx++=dc*w+ds*z;
			*dy++=-ds*w+dc*z;
		}
	}
}
