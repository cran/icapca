/* Copyright 2012-2014 Roger P. Woods, M.D. */
/* Modified 6/11/2014 */

/*
 * compute Y=R*X where X and Y are N X M matrices and R is a random orthonormal matrix
 *
 */
 
#include "AIR.h"
#ifdef USING_R
#include <R.h>
#include <Rmath.h>
#endif


static void mulmat(double **a, const unsigned int m, const unsigned int n, double **b, const unsigned int r, double **c)
{
	/* a is mxn */
	/* b is nxr */
	/* c is mxr */
	unsigned int y;
	for(y=0;y<r;y++){
		unsigned int x;
		
		for(x=0;x<m;x++){
			c[y][x]=0.0;
			{
				unsigned int i;
				for(i=0;i<n;i++){
					c[y][x]+=a[i][x]*b[y][i];
				}
			}
		}
	}
}

void AIR_random_rotation(const unsigned int skip, double **r, const unsigned int n, const unsigned int m, double **x, double **y, double ***wrk)
{
	{
		double pi=2.0*acos(0);
		unsigned int randomization;
		
#ifdef USING_R
		GetRNGstate();
#endif
		
		for(randomization=0; randomization<=skip; randomization++){
			AIR_Error errcode=1;
			do{
				unsigned int i;
				for(i=0; i<n; i++){
	
					r[i][i]=0.0;
	
					{
						unsigned int j;
						for(j=0; j<i; j++){
						
#ifndef USING_R
							double angle=2*pi*rand()/(RAND_MAX+1.0);
#else
							double angle=2*pi*unif_rand();
#endif
							if(angle>pi) angle-=2*pi;
							
							r[i][j]=angle;
							r[j][i]=-angle;
						}
					}
				}
				errcode=AIR_eexper_pade(n, r, wrk, 0);
			} while (errcode!=0);
		}
#ifdef USING_R
		PutRNGstate();
#endif
	}
	mulmat(r, n, n, x, m, y);
}
