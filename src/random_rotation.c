/* Copyright 2012 Roger P. Woods, M.D. */
/* Modified 10/12/2012 */

/*
 * compute Y=R*X where X and Y are N X M matrices and R is a random orthonormal matrix
 *
 */
 
 #include "AIR.h"

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
		
		for(randomization=0; randomization<=skip; randomization++){
			AIR_Error errcode=1;
			do{
				unsigned int i;
				for(i=0; i<n; i++){
	
					r[i][i]=0.0;
	
					{
						unsigned int j;
						for(j=0; j<i; j++){
						
							double angle=2*pi*rand()/(RAND_MAX+1.0);
							if(angle>pi) angle-=2*pi;
							
							r[i][j]=angle;
							r[j][i]=-angle;
						}
					}
				}
				errcode=AIR_eexper_pade(n, r, wrk, 0);
			} while (errcode!=0);
		}
	}
	mulmat(r, n, n, x, m, y);
}
