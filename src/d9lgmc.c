/* Copyright 1995-2013 Roger P. Woods, M.D. */
/* Modified 7/4/13 */

/* double d9lgmc()
 *	
 * computes the log gamma correction factor for X>=10, so that
 * log(dgamma(x))=log(sqrt(2*pi))+(x-.5)*log(x)-x+d9lgmc(x)		
 */


#include "AIR.h"
#include <float.h>

double AIR_d9lgmc(const double x, AIR_Error *errcode)

{
	static AIR_Boolean first=TRUE;
	static double algmcs[15]={.1666389480451863247205729650822e+0,-.1384948176067563840732986059135e-4,.9810825646924729426157171547487e-8,-.1809129475572494194263306266719e-10,.6221098041892605227126015543416e-13,-.3399615005417721944303330599666e-15,.2683181998482698748957538846666e-17,-.2868042435334643284144622399999e-19,.3962837061046434803679306666666e-21,-.6831888753985766870111999999999e-23,.1429227355942498147573333333333e-24,-.3547598158101070547199999999999e-26,.1025680058010470912000000000000e-27,-.3401102254316748799999999999999e-29,.1276642195630062933333333333333e-30};
	static double xbig;
	static double xmax;
	static unsigned int nalgm;

	if(first){
		nalgm=AIR_initds(algmcs,15,DBL_EPSILON/FLT_RADIX);
		if(nalgm==0){
#ifndef USING_R
			printf("%s: %d: ",__FILE__,__LINE__);
			printf("d9lgmc error initializing\n");
#else
			REprintf("%s: %d: ",__FILE__,__LINE__);
			REprintf("d9lgmc error initializing\n");
#endif
			*errcode=AIR_D9LGMC_CANT_INIT_ERROR;
			return 0.0;
		}
		xbig=1.0/sqrt(DBL_EPSILON/FLT_RADIX);
		{
			double temp=log(DBL_MAX/12.0);
			{
				double temp2=-log(12.0*DBL_MIN);
				if(temp>temp2) temp=temp2;
			}
			xmax=exp(temp);
		}
		first=FALSE;
	}

	if(x<10.0){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("d9lgmc error, x must be>=10\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("d9lgmc error, x must be>=10\n");
#endif
		*errcode=AIR_D9LGMC_SMALL_X_ERROR;
		return 0.0;
	}
	if(x>=xmax){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("WARNING: d9lgmc, x so big that value underflows\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("WARNING: d9lgmc, x so big that value underflows\n");
#endif
		return 0.0;
	}
	if(x<xbig){
		double answer=AIR_dcsevl(2.0*(10.0/x)*(10.0/x)-1.0,algmcs,nalgm,errcode);
		return answer/x;
	}
	return 1.0/(12.0*x);
}

	
