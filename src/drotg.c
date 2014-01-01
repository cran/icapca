/* Copyright 1996-99 Roger P. Woods, M.D. */
/* Modified 11/19/99 */

/*
 * void drotg()
 */


#include "AIR.h"

void AIR_drotg(double *da, double *db, double *dc, double *ds)

{
	if(fabs(*da)>fabs(*db)){

		double u=*da+*da;
		double v=*db/u;
		double r=sqrt(.25+v*v)*u;

		*dc=*da/r;
		*ds=v*(*dc+*dc);
		*db=*ds;
		*da=r;
		return;
	}
	if(*db==0.0){
		*dc=1.0;
		*ds=0.0;
		return;
	}
	{
		double u=*db+*db;
		double v=*da/u;

		*da=sqrt(.25+v*v)*u;
		*ds=*db/(*da);
		*dc=v*(*ds+*ds);
	}
	if(*dc!=0.0){
		*db=1.0/(*dc);
		return;
	}
	*db=1.0;
	return;
}
