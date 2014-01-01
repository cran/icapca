/* Copyright 2003-2012 Roger P. Woods, M.D. */
/* Modified 10/9/12 */
      
/*
 * WARNING: the array wrk must be allocated using wrk=AIR_matrix3(m,m,>=4) with m matching
 * the value of m used to call this routine.
 *
 * Compute the exponential of the matrix H using the scaling and squaring method
 * with Pade approximation. The method is implemented as described on pages 134-136 of
 * Sidje, R. Expokit: A software package for computing matrix exponentials. ACM
 * Transactions on Mathematical Software 1998;24(1):130-156.
 *
 * On Input:
 *
 * 	M is the order of the input matrix H
 * 	H is the input matrix
 * 	WRK is a work array allocated using matrix3(m,m,>=4);
 *
 * Returns:
 * 	O if no errors
 * 	errcode otherwise
 *
 * On Output:
 *
 * 	H contains the matrix exponential of the input H
 *
 */

#include "AIR.h"
#define DEGREE 6	// Degree of Pade approximation
//#define MAXIMUM_MATRIX_DIMENSION 10

AIR_Error AIR_eexper_pade(const unsigned int m, double **h, double ***wrk, AIR_Boolean is_affine)
{
	//unsigned int mm=m*m;
	double coefficients[DEGREE+1];
	
	double **h2=wrk[0];
	double **p=wrk[1];
	double **q=wrk[2];
	double **work=wrk[3];
	
	unsigned int ns;
	unsigned int odd;

	// Find ns such that ||H/2^ns|| < 1/2
	
	{
		unsigned int i;
		for(i=0;i<m;i++){
			wrk[0][0][i]=0.0;
		}
	}
	{
		unsigned int j;
		for(j=0;j<m;j++){
			unsigned int i;
			for(i=0;i<m;i++){
				wrk[0][0][i]+=fabs(h[j][i]);
			}
		}
	}
	{
		double norm=0.0;
		{
			unsigned int i;
			for(i=0;i<m;i++){
				if(wrk[0][0][i]>norm) norm=wrk[0][0][i];
			}
		}
		if(norm==0.0){
			unsigned int j;
			for(j=0;j<m;j++){
				h[j][j]=1.0;
			}
			return 0;
		}
		{
			double temp=log(norm);
			temp/=log(2.0);
			{
				int ns1=(int)temp+2;
				if(0>ns1) ns1=0;
				ns=(unsigned int)ns1;
			}
		}
	}
	{
		double scale=1.0/pow(2.0,(double)ns);
		double scale2=scale*scale;
	
		// Compute Pade coefficients
		{
			unsigned int i=DEGREE+1;
			unsigned int j=2*DEGREE+1;
			coefficients[0]=1.0;
			{
				unsigned int k;
				for(k=1;k<=DEGREE;k++){
					coefficients[k]=(coefficients[k-1]*(double)(i-k))/(double)(k*(j-k));
				}
			}
		}
	
		// H2=scale2*H*H
		
		AIR_dgemm('n','n',m,m,m,scale2,h,h,0.0,h2);
	
		// Initialize p (numerator) and q (denominator)
		{
			double cp=coefficients[DEGREE-1];
			double cq=coefficients[DEGREE];
			{
				unsigned int j;
				for(j=0;j<m;j++){
					unsigned int i;
					for(i=0;i<m;i++){
						p[j][i]=0.0;
						q[j][i]=0.0;
					}
					p[j][j]=cp;
					q[j][j]=cq;
				}
			}
		}
		
		// Horner rule
		odd=1;
		{
			unsigned int k=DEGREE-1;
			do{
				double **current;
				if(odd!=0) current=q;
				else current=p;
				
				AIR_dgemm('n','n',m,m,m,1.0,current,h2,0.0,work);
				{
					unsigned int j;
					
					for(j=0;j<m;j++){
						work[j][j]+=coefficients[k-1];
					}
				}
				if(odd!=0) q=work;
				else p=work;
				
				work=current;
				odd=1-odd;
				k--;
			}while(k>0);
		}
		
		
		// (+/-)(1+2*(p\q))
		if(odd==1){
			AIR_dgemm('n','n',m,m,m,scale,q,h,0.0,work);
			q=work;
		}
		else{
			AIR_dgemm('n','n',m,m,m,scale,p,h,0.0,work);
			p=work;
		}
	}
	{
		unsigned int j;
		for(j=0;j<m;j++){
			AIR_daxpy(m,-1.0,p[j],q[j]);
		}
	}
	//AIR_daxpy(mm,-1.0,p[0],q[0]);
	
	// Solve inv(q)*p 
	{
		unsigned int ipvt[m]; // Previously dimensioned as MAXIMUM_MATRIX_DIMENSION
		
		if(AIR_dgefa(q,m,ipvt)!=m){
			return AIR_SINGULAR_EEXPER_ERROR;
		}
		else{
			unsigned int j;
			
			for(j=0;j<m;j++){
				AIR_dgesl(q,m,ipvt,p[j],0);
			}
		}
	}
	{
		unsigned int j;
		for(j=0;j<m;j++){
			AIR_dscal(m,2.0,p[j]);
		}
	}
	//AIR_dscal(mm,2.0,p[0]);
	{
		unsigned int j;
		for(j=0;j<m;j++){
			p[j][j]+=1.0;
		}
	}
	{
		double **out=p;
		
		if(ns==0 && odd==1){
			{
				unsigned int j;
				for(j=0;j<m;j++){
					AIR_dscal(m,-1.0,p[j]);
				}
			}
			//AIR_dscal(mm, -1.0, p[0]);
		}
		else{
	
			// Squaring
			odd=1;
			{
				unsigned int k;
				
				for(k=0;k<ns;k++){
				
					double **in;
					
					if(odd!=0){
						in=p;
						out=q;
					}
					else{
						in=q;
						out=p;
					}
					AIR_dgemm('n','n',m,m,m,1.0,in,in,0.0,out);
					odd=1-odd;
				}
			}
		}
		{
			unsigned int j;
			for(j=0;j<m;j++){
				unsigned int i;
				for(i=0;i<m;i++){
					h[j][i]=out[j][i];
				}
			}
		}
	}
	if(is_affine) h[m-1][m-1]=1.0;
	
	return 0;
}
