/* Copyright 2004-2010 Roger P. Woods, M.D. */
/* Modified 10/08/10 */

/* void dqrdc()
 *
 * Uses Householder transformations to compute the qr factorization of an n by p matrix x.
 * Column pivoting based on the 2-norms of the reduced columns may be performed at the users option.
 *
 *	On Entry
 *
 *	x		double[p][n]	matrix to decompose
 *	n 		integer			number of rows of x
 *	p		integer			number of columes of x
 *	jpvt	integer[p]		integers that control the selection of pivot columns
 *								if jpvt[k] >0 x[k] is an initial column
 *								if jpvt[k] ==0 x[k] is a free column
 *								if jpvt[k] <0, x[k] is a final column
 *
 *
 *							Pivoting is useful when the matrix has less than full rank because
 *							it allows vectors that span the full subspace to be characterized in the
 *							initial columns while columns that are just linear combinations of these
 *							are moved to the latter columns. For these latter trivial components,
 *							the corresponding diagonal element of R becomes negligible
 *
 *							Note that marking a column as an initial column allows that column
 *							to be forced into a least squares file.
 *
 *							Marking a column as a final column is appropriate when working with the
 *							augmented matrix (X with a tilde over it).
 *							
 *							jpvt is only used if job!=0
 *
 *	work	double[p]		work array used only if job!=0
 *	job		integer
 *				if job==0, no pivoting (see explanation of jpvt above)
 *				if job!=0, pivoting is done
 *
 *	On Exit
 *	
 *	x		contains in its upper triangle the upper triangular matrix r of the qr 
 *			factorization. Below its diagonal, it contains information from which the orthogonal
 *			part of the decomposition can be recovered. Note that if pivoting has been requested,
 *			the decomposition is not that of the original matrix x, but that of x with its columns
 *			permuted as described by jpvt
 *
 *	qraux	double[p] countains further information required to recover the orthogonal part of the
 *			decomposition
 *
 *	jpvt	jpvt[k] contains the index of the column of the original matrix that has been 
 *			interchanged into the k-th column, if pivoting was requested.
 *
 *	Use dqrsl to recover Q by computing Qy for properly chosen values of y
 */
 
 /* Note that values in jvpt follow Fortran convention until immediately before the
  *  function returns. Otherwise, it would not be possible to mark the first column
  *  as belonging at the end.
  */
 
#include "AIR.h"

void AIR_dqrdc(double **x, const unsigned int n, const unsigned int p, int *jpvt, double *work, const unsigned int job, double *qraux)
{
	int pl=0;
	int pu=-1;
	if(job!=0){
		// pivoting requested, rearrange columns according to jpvt
		{
			int j;
			for(j=0;j<p;j++){
				int swapj=jpvt[j]>0;
				int negj=jpvt[j]<0;
				jpvt[j]=j+1;
				if(negj) jpvt[j]=-(j+1);
				if(swapj){
					if(j!=pl) AIR_dswap(n,x[pl],x[j]);
					jpvt[j]=jpvt[pl];
					jpvt[pl]=j+1;
					pl++;
				}
			}
		}
		pu=p-1;
		{
			int j;
			for(j=p-1;j>=0;j--){
				if(jpvt[j]<0){
					jpvt[j]=-jpvt[j];
					if(j!=pu){
						AIR_dswap(n,x[pu],x[j]);
						{
							int jp=jpvt[pu];
							jpvt[pu]=jpvt[j];
							jpvt[j]=jp;
						}
					}
					pu--;
				}
			}
		}
	}

	// Compute the norms of the free columns
	if(pu>=pl){
		int j;
		for(j=pl;j<=pu;j++){
			qraux[j]=AIR_dnrm2(n,x[j]);
			work[j]=qraux[j];
		}
	}
	// Perform the Householder reduction of x
	{	
		int lup=n;
		int l;
				
		if(p<n) lup=p;

		for(l=0;l<lup;l++){

			if(l>=pl && l<pu){

				// Locate the column of largest norm and bring it into the pivot position
				double maxnrm=0.0;
				int maxj=l;
				{
					{
						int j;
						for(j=l;j<=pu;j++){
							if(qraux[j]>maxnrm){
								maxnrm=qraux[j];
								maxj=j;
							}
						}
					}
					if(maxj!=l){
						AIR_dswap(n,x[l],x[maxj]);
						qraux[maxj]=qraux[l];
						work[maxj]=work[l];
						{
							int jp=jpvt[maxj];
							jpvt[maxj]=jpvt[l];
							jpvt[l]=jp;
						}
					}
				}
			}
			qraux[l]=0.0;
			if(l!=n-1){

				// Compute the Householder transformation for column l

				double nrmxl=AIR_dnrm2(n-l,&x[l][l]);

				if(nrmxl!=0.0){
				
					if(x[l][l]!=0.0){
						if(x[l][l]>0.0){
							nrmxl=fabs(nrmxl);
						}
						else{
							nrmxl=-fabs(nrmxl);
						}
					}
					AIR_dscal(n-l,1.0/nrmxl,&x[l][l]);

					x[l][l]+=1.0;
		
					// Apply the transformation to the remaining columns, updating the norms
					{
						int lp1=l+1;
						if(p-1>=lp1){
							int j;
	
							for(j=lp1;j<p;j++){
								double t=-AIR_ddot(n-l,&x[l][l],&x[j][l])/x[l][l];
	
								AIR_daxpy(n-l,t,&x[l][l],&x[j][l]);
	
								if(j>=pl && j<=pu){
								
									if(qraux[j]!=0.0){
										double temp=fabs(x[j][l])/qraux[j];
										double tt=1.0-temp*temp;
										
										if(0.0>tt) tt=0.0;
										t=tt;
										temp=qraux[j]/work[j];
										tt=1.0+0.05*tt*temp*temp;
	
										if(tt!=1.0){
											qraux[j]=qraux[j]*sqrt(t);
										}
										else{
											qraux[j]=AIR_dnrm2(n-l-1,&x[j][l+1]);
											work[j]=qraux[j];
										}
									}
								}
							}
						}
					}

					// Save the transformation
					qraux[l]=x[l][l];
					x[l][l]=-nrmxl;
				}
			}
		}
	}
	if(job!=0){
		unsigned int i;
		for(i=0; i<p; i++){
			jpvt[i]--;
		}
	}
}
