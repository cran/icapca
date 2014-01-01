/* Copyright 2004-2005 Roger P. Woods, M.D. */
/* Modified 2/15/05 */

/*	void dqrsl()
 *
 *	On entry:
 *		x		double[p][n]	the output of dqrdc
 *		n		number of rows in matrix xk. Must be the same as n in dqrdc
 *		k		number of columns in matrix xk. It must not be greater than min(n,p)
 *				where p is the same as in calling dqrdc
 *		draux	double[p]	from dqrdc
 *		y		double[n] n-vector to be manipulated by dqrsl
 *		job		int	
 *					Has the decimal expansion abcde:
 *						if a!=0, compute qy
 *						if b, c, d or e !=0, compute qty
 *						if c!=0, compute b
 *						if d!=0, compute rsd
 *						if e!=0, compute xb
 *					
 *					note that a request to compute b, rsd or xb automatically
 *					triggers the computation of qty, for which an array must be
 *					provided in the calling sequence
 *
 *	On return:
 *		qy		double[n]		contains q*y, if requested
 *		qty		double[n]		contains trans(q)*y, if its computation has been requested
 *		b		double[k]		contains the solution of the least squares problem:
 *									minimize norm2(y-xk*b), if its computation has been requested
 *									(note that if pivoting was requested in dqrdc, the j-th
 *									component of b will be associated with column jpvt[j]
 *									of the original matrix x that was input into dqrdc
 *								
 *		rsd		double[n]		contains the least squares residual y-xk*b, if requested
 *								It is also the orthogonal projection of y onto the orthogonal
 *								complement of the column space of xk
 *							
 *		xb		double[n]		xb contains the least squares approximation xk*b, if requested
 *								It is also the orthogonal projection of y on the column space of x
 *						
 *		info	int				-1, unless computation of b has been requested and r is exactly
 *								singular. In this case, info is the index of the first zero
 *								diagonal element of r and b is left unaltered
 */
 
#include "AIR.h"

static void dcopy(unsigned int n, double *a, double *b){

	unsigned int i;
	for(i=0;i<n;i++){
		*b++=*a++;
	}
}

void AIR_dqrsl(double **x, unsigned int n, unsigned int k, double *qraux, double *y, int job, double *qy, double *qty, double *b, double *rsd, double *xb, int *info)
{
	*info=-1;
	{
		// Determine what  is to be computed
		int cqy=job/10000!=0;
		int cqty=job%10000!=0;
		int cb=(job%1000)/100!=0;
		int cr=(job%100)/10!=0;
		int cxb=(job%10)!=0;
		int ju=n-1;
		if(k<n-1) ju=k;
		
		if(ju==0){
			// Special action when n=1
			if(cqy) qy[0]=y[0];
			if(cqty) qty[0]=y[0];
			if(cxb) xb[0]=y[0];
			if(cb){
				if(x[0][0]==0.0){
					*info=0;
				}
				else{
					b[0]=y[0]/x[0][0];
				}
			}
			if(cr) rsd[0]=0.0;
		}
		else{
			// Set up to compute qy or qty
			if(cqy) dcopy(n,y,qy);
			if(cqty) dcopy(n,y,qty);
			if(cqy){
				// Compute qy
				int j;
				for(j=ju-1;j>=0;j--){
					if(qraux[j]!=0.0){
						double temp=x[j][j];
						x[j][j]=qraux[j];
						{
							double t=-AIR_ddot(n-j,&x[j][j],&qy[j])/x[j][j];
							AIR_daxpy(n-j,t,&x[j][j],&qy[j]);
						}
						x[j][j]=temp;
					}
				}
			}
			if(cqty){
				// Compute trans(q)*y
				int j;
				for(j=0;j<ju;j++){
					if(qraux[j]!=0.0){
						double temp=x[j][j];
						x[j][j]=qraux[j];
						{
							double t=-AIR_ddot(n-j,&x[j][j],&qty[j])/x[j][j];
							AIR_daxpy(n-j,t,&x[j][j],&qty[j]);
						}
						x[j][j]=temp;
					}
				}
			}
			// Set up to compute b, rsd or xb
			if(cb) dcopy(k,qty,b);
			{
				int kp1=k+1;
				if(cxb) dcopy(k,qty,xb);
				if(cr && (k<n)) dcopy(n-k, &qty[k], &rsd[k]);
				if(cxb && (kp1<=n)){
					unsigned int i;
					for(i=k;i<n;i++){
						xb[i]=0.0;
					}
				}
				if(cr){
					unsigned int i;
					for(i=0;i<k;i++){
						rsd[i]=0.0;
					}
				}
				if(cb){
					// compute b
					int j;
					for(j=k-1;j>=0;j--){
						if(x[j][j]==0.0){
							*info=j;
							break;
						}
						b[j]=b[j]/x[j][j];
						if(j!=0){
							double t=-b[j];
							AIR_daxpy(j,t,x[j],b);
						}
					}
				}
				if(cr || cxb){
					// Compute rsd or xb as required
					int j;
					for(j=ju-1;j>=0;j--){
						if(qraux[j]!=0.0){
							double temp=x[j][j];
							x[j][j]=qraux[j];
							if(cr){
								double t=-AIR_ddot(n-j,&x[j][j],&rsd[j])/x[j][j];
								AIR_daxpy(n-j,t,&x[j][j],&rsd[j]);
							}
							if(cxb){
								double t=-AIR_ddot(n-j,&x[j][j],&xb[j])/x[j][j];
								AIR_daxpy(n-j,t,&x[j][j],&xb[j]);
							}
							x[j][j]=temp;
						}
					}
				}
			}
		}
	}
}
