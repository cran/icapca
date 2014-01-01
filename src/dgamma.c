/* Copyright 1995-2013 Roger P. Woods, M.D. */
/* Modified 7/4/13 */

/*
 * double dgamma()
 *
 * calculates the complete gamma function for X
 */

#include "AIR.h"
#include <float.h>

double AIR_dgamma(const double x, AIR_Error *errcode)

{
	static const double gamcs[42]={.8571195590989331421920062399942e-2,.4415381324841006757191315771652e-2,.5685043681599363378632664588789e-1,-.4219835396418560501012500186624e-2,.1326808181212460220584006796352e-2,-.1893024529798880432523947023886e-3,.3606925327441245256578082217225e-4,-.6056761904460864218485548290365e-5,.1055829546302283344731823509093e-5,-.1811967365542384048291855891166e-6,.3117724964715322277790254593169e-7,-.5354219639019687140874081024347e-8,.9193275519859588946887786825940e-9,-.1577941280288339761767423273953e-9,.2707980622934954543266540433089e-10,-.4646818653825730144081661058933e-11,.7973350192007419656460767175359e-12,-.1368078209830916025799499172309e-12,.2347319486563800657233471771688e-13,-.4027432614949066932766570534699e-14,.6910051747372100912138336975257e-15,-.1185584500221992907052387126192e-15,.2034148542496373955201026051932e-16,-.3490054341717405849274012949108e-17,.5987993856485305567135051066026e-18,-.1027378057872228074490069778431e-18,.1762702816060529824942759660748e-19,-.3024320653735306260958772112042e-20,.5188914660218397839717833550506e-21,-.8902770842456576692449251601066e-22,.1527474068493342602274596891306e-22,-.2620731256187362900257328332799e-23,.4496464047830538670331046570666e-24,-.7714712731336877911703901525333e-25,.1323635453126044036486572714666e-25,-.2270999412942928816702313813333e-26,.3896418998003991449320816639999e-27,-.6685198115125953327792127999999e-28,.1146998663140024384347613866666e-28,-.1967938586345134677295103999999e-29,.3376448816585338090334890666666e-30,-.5793070335782135784625493333333e-31};
	static double pi;
	static double sq2pil;
	static double 	dxrel;
	static unsigned int ngam;
	static double 	xmin,xmax;
	static AIR_Boolean first=(AIR_Boolean)1;	

	*errcode=0;

	if(first){
		pi=2.0*acos(0.0);
		sq2pil=log(sqrt(2.0*pi));
		dxrel=sqrt(DBL_EPSILON);
		ngam=AIR_initds(gamcs,42,0.1*DBL_EPSILON/FLT_RADIX);
		
		*errcode=AIR_dgamlm(&xmin,&xmax);
		if(*errcode!=0) return 0.0;
		
		first=0;
	}
	{
		double y=fabs(x);
		
		if(y<=10){

			/*Compute gamma(x) for -xbnd<=x<=xbnd*/
			/*Reduce interval and find gamma(1+y) for */
			/* 0.0<=y<1.0 first of all*/
			{
				int n=(int)x;
				
				if(x<0.0) n--;
				y=x-n;
				n--;
				{
					double answer=AIR_dcsevl(2.0*y-1.0,gamcs,ngam,errcode);
					if(*errcode!=0) return 0.0;
										
					answer+=.9375;
					if(n==0) return answer;
					
					if(n<=0){
						/* Compute gamma(x) for x<1.0 */
						n=-n;
						if(x==0.0){
#ifndef USING_R
							printf("%s: %d: ",__FILE__,__LINE__);
							printf("dgamma: x is 0\n");
#else
							REprintf("%s: %d: ",__FILE__,__LINE__);
							REprintf("dgamma: x is 0\n");
#endif
							*errcode=AIR_DGAMMA_X_IS_ZERO_ERROR;
							return 0.0;
						}
						if(x<0.0 && x+n-2==0.0){
#ifndef USING_R
							printf("%s: %d: ",__FILE__,__LINE__);
							printf("dgamma: x is negative integer\n");
#else
							REprintf("%s: %d: ",__FILE__,__LINE__);
							REprintf("dgamma: x is negative integer\n");
#endif
							*errcode=AIR_GAMMA_OF_NEGATIVE_ERROR;
							return 0.0;
						}
						if(x<-.5 && fabs((x-ceil(x-.5))/x)<dxrel){
#ifndef USING_R
							printf("%s: %d: ",__FILE__,__LINE__);
							printf("dgamma: WARNING answer < half precision because x too near negative integer\n");
#else
							REprintf("%s: %d: ",__FILE__,__LINE__);
							REprintf("dgamma: WARNING answer < half precision because x too near negative integer\n");
#endif
						}
						{
							int i;
							
							for(i=1;i<=n;i++){
							
								answer/=(x+i-1);

							}
						}
					}
					else{
						/* gamma(x) for x>=2.0 and x <=10.0 */

						int i;
						
						for(i=1;i<=n;i++){
							answer*=(y+i);
						}
					}

					return answer;
				}
			}
		}
		else{		
			/* gamma(x) for fabs(x) >10.0. Recall y=fabs(x) */
			if(x>xmax){
#ifndef USING_R
				printf("%s: %d: ",__FILE__,__LINE__);
				printf("dgamma: x so big gamma overflows\n");
#else
				REprintf("%s: %d: ",__FILE__,__LINE__);
				REprintf("dgamma: x so big gamma overflows\n");
#endif
				*errcode=AIR_GAMMA_OVERFLOW_ERROR;
				return 0.0;
			}
			
			if(x<xmin){
#ifndef USING_R
				printf("%s: %d: ",__FILE__,__LINE__);
				printf("dgamma: WARNING x so small gamma underflows\n");
#else
				REprintf("%s: %d: ",__FILE__,__LINE__);
				REprintf("dgamma: WARNING x so small gamma underflows\n");
#endif
				return 0.0;
			}
			{
				double answer=exp((y-.5)*log(y)-y+sq2pil+AIR_d9lgmc(y,errcode));
				if(*errcode!=0) return 0.0;

				if(x>0) return answer;
				
				if(fabs((x-ceil(x-.5))/x)<dxrel){
#ifndef USING_R
					printf("%s: %d: ",__FILE__,__LINE__);
					printf("dgamma: WARNING answer < half precision, x too near negative integer\n");
#else
					REprintf("%s: %d: ",__FILE__,__LINE__);
					REprintf("dgamma: WARNING answer < half precision, x too near negative integer\n");
#endif
				}
				{
					double sinpiy=sin(pi*y);
					
					if(sinpiy==0.0){
#ifndef USING_R
						printf("%s: %d: ",__FILE__,__LINE__);
						printf("dgamma: x is a negative integer\n");
#else
						REprintf("%s: %d: ",__FILE__,__LINE__);
						REprintf("dgamma: x is a negative integer\n");
#endif
						*errcode=AIR_GAMMA_OF_NEGATIVE_ERROR;
						return 0.0;
					}
					return -pi/(y*sinpiy*answer);
				}
			}
		}
	}
}

	
