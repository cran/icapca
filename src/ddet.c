/* Copyright 1996-2013 Roger P. Woods, M.D. */
/* Modified 02/01/13 */


/* 
 * double ddet()
 *
 * Compute the determinant of matrix of order n using factors computed by
 *  DGECO or DGEFA
 * 
 * This is adapted from DGEDI
 *
 * On entry
 *	a	the output from DGECO or DGEFA
 *	n	the order of the matrix A
 *	*ipvt	the pivot vector from DGECO or DGEFA
 *
 * Returns
 *	determinant of original matrix
 */

#include "AIR.h"

double AIR_ddet(double **a, const unsigned int n, unsigned int *ipvt)

{
	double	det[2]={1.0, 0.0};
	unsigned int i;

	for(i=0;i<n;i++){
		if(ipvt[i]!=i) det[0]*=-1.0;
		det[0]*=a[i][i];
		if(det[0]==0.0) break;
		while(fabs(det[0])<1.0){
			det[0]*=10.0;
			det[1]-=1.0;
		}
		while(fabs(det[0])>=10.0){
			det[0]/=10.0;
			det[1]+=1.0;
		}
	}
	return (det[0]*pow(10.0,det[1]));
}
