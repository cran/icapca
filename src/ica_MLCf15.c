/* Copyright 2004-2013 Roger P. Woods, M.D. */
/* Modified: 7/4/13 */

/*
 * subroutine for ICA minimization with sub-Gaussian, super-Gaussian and Gaussian sources
 *
 * Generalized subroutine for ICA computation with:
 *	Sub-Gaussian sources
 *	Super-Gaussian sources
 *	Gaussian soures orthogonal to all non-Gaussian sources
 *
 *	IMPORTANT: This routine has two, mutually incompatible operative modes:
 *		Mode 1: dw_temp==NULL	returned sources are valid, gradients are NOT computed
 *		Mode 2: dw_temp!=NULL	returned sources are NOT valid, gradients are computed
 *	
 *
 * Returns the negative log likelihood and its gradient with respect to W
 *
 * On input,
 *	UNUSED	an integer that is ignored (included for compatibility with ucminf)
 *	W_PARAMETERS is the set of parameters being optimized (not changed by this routine)
 *  dW_TEMP is space for storing the gradients, if NULL, derivatives won't be computed
 *	dW is the set of computed gradients of W
 *	CONSTANTS is a set of parameter pointers:
 *		[0] 	*M is the number of rows in the matrix to be factored
 *		[1] 	*N is the number of columns in the matrix to be factored
 *		[2] 	X is the M by N matrix to be factored into independent components
 *		[3] 	W is an M by M matrix 
 *		[4] 	S is the M by N matrix of computed sources (valid only if dw_temp==NULL)
 *		[5] 	IPVT is an unsigned int array of length M
 *		[6] 	dW is an M by M double array
 *		[7] 	WORK is a pointer to a double array of length M
 *		[9] 	SUBGC is the number of sub-Gaussian components
 *		[10] 	SUPERGC is the number super-Gaussian components
 *		[11] 	GC is the number of Gaussian components
 *		[14] 	QR_MATRIX is an M by M double matrix
 *		[15] 	QRAUX and E is a double vector of length M
 *		[16] 	NULL_MATRIX_COLUMN and D is a double vector of length N  ****previously length M****
 *		[17] 	SVD_WORK is a double vector of length N
 *		[18] 	S2 is an N by M matrix
 *		[19]	VARIANCE is a double vector of length M
 *
 *	On output,
 *		the negative log likelihood is returned
 *		dW contains pointers to dW_TEMP
 *		if not NULL, dW_TEMP contains the derivatives of the negative log likelihood with respect to W
 *		
 */

#include "AIR.h"

static double logcosh(double x){
	double result=log(cosh(x));
	if(isinf(result)) result=fabs(x)-log(2);
	return(result);
}

double AIR_ica_MLCf15(/*@unused@*/ const unsigned int unused, const double *w_parameters, double *dw_temp, void **constants){
	
	static AIR_Boolean initialized=FALSE;
	static double pi;
	static double log_pi;
	static double log_2pi;
	static double log_2;
	static double sqrt_2;
	
	const unsigned int m=*((unsigned int *)constants[0]);	/* m<=n */
	const unsigned int n=*((unsigned int *)constants[1]);	/* m<=n */
	
	double **x=constants[2];
	double **w=constants[3];
	double **s=constants[4];
	
	const unsigned int subgc=*((unsigned int *)constants[9]);			/* Number of sub-gaussian components */
	const unsigned int supergc=*((unsigned int *)constants[10]);		/* Number of super-gaussian components */
	const unsigned int gc=m-subgc-supergc;				/* Number of Gaussian components */

	const unsigned int ngc=subgc+supergc; /* Number of non-gaussian components */
	
	double *logdeterminant=(double *)constants[11];		/* log of absolute value of determinant is stored here */
	
	double f=0.0;
	
#ifdef USING_R
	R_CheckUserInterrupt();
#endif
#ifdef USING_MATLAB
	/* This doesn't work on all platforms, but may allow
	 * Ctrl-C to interrupt execution on some platforms
	 */
	mexEvalString("drawnow;");
#endif
		
	if(!initialized){
		pi=2.0*acos(0.0);
		log_pi=log(pi);
		log_2pi=log(2.0*pi);
		log_2=log(2.0);
		sqrt_2=sqrt(2.0);
		initialized=TRUE;
	}
	/*
	 * Compute W, the ICA unmixing matrix
	 *
	 * The next ncg rows of w, corresponding to subgaussian and
	 *  supergaussian sources, are extracted from w_parameters
	 *
	 * The remaining rows of w are constructed as an orthonormal 
	 *  basis of the subspace that is orthogonal to the
	 *  non-gaussian sources
	 */
	{
		double **qr_matrix=constants[14];	/* m by m  */

		/* Copy adjustable non-gaussian parameters from w_parameters to w
		 *
		 *  w(1:ngc,:)=w_parameters(1:ngc,:)
		 *
		 * and transpose these elements of w into qr_matrix
		 *
		 *  qr_matrix(:,1:ncg)=w(1:ngc,:)'
		 */
		{
			const double *temp=w_parameters;
			unsigned int j;
					
			for(j=0; j<m; j++){
			
				unsigned int i;
				
				for(i=0 ; i<ngc; i++){
				
					qr_matrix[i][j]=*temp;
					w[j][i]=*temp++;
				}
			}
		}
		if(ngc<m){
			{
				/* Perform QR of qr_matrix, replacing qr_matrix with R
				 *
				 * [q, r]=qr(qr_matrix(:,1:ngc))
				 * qr_matrix(:,1:ngc)=r
				 */
				double *qraux=constants[15];		/* vector of length ngc (actually length m) */

				AIR_dqrdc(qr_matrix, m, ngc, NULL, NULL, 0, qraux);
				{
					/* Compute remaining rows of w to span orthogonal subspace
					 *
					 *  w(ngc+1:end,:)=q(:,ngc+1:end)'
					 *
					 * Note that w is set to an identity matrix if ngc==0
					 */
				 
					double *null_matrix_column=constants[16];	/* vector of length m (actually of length n) */
					unsigned int i;
					
					for(i=ngc; i<m; i++){
						{
							unsigned int j;
							
							for(j=0; j<m; j++){
								if(i==j) null_matrix_column[j]=1.0;
								else null_matrix_column[j]=0.0;
							}
						}
						{
							int info;	/* No need to check afterwards when job==10000 */
							int job=10000;
							AIR_dqrsl(qr_matrix, m, ngc, qraux, null_matrix_column, job, null_matrix_column, NULL, NULL, NULL, NULL, &info);
						}
						{
							unsigned int j;
							
							for(j=0; j<m; j++){
								w[j][i]=null_matrix_column[j];
							}
						}
					}
				}
			}
			/* Compute the non-gaussian sources and the subspace
			 *  holders for the gaussian sources
			 *
			 * s=w*x
			 */
			AIR_dgemm('n', 'n', m, n, m, 1.0, w, x, 0.0, s);
			{
				double **s2=constants[18];			/* n by m */
				{
					/* Transpose s into s2
					 *
					 * s2=s'
					 */
					unsigned int j;
					
					for(j=0; j<n; j++){
					
						unsigned int i;
						
						for(i=0; i<m; i++){
						
							s2[i][j]=s[j][i];
						}
					}
				}
				{
					/* Perform QR decomposition of s2
					 *
					 * [q, r]=qr(s2,0)
					 */
					double *qraux=constants[15]; 	/* vector of length m */

					AIR_dqrdc(s2, n, m, NULL, NULL, 0, qraux);
					{
						/* Compute gaussian sources
						 *
						 * s(ngc+1:end,:)=sqrt(n)*q(:,ngc+1:end)'
						 */
						double *null_matrix_column=constants[16];	/* vector of length n */
						double sqrt_n=sqrt(n);

						unsigned int i;

						for(i=ngc; i<m; i++){
							{
								unsigned int j;
								for(j=0; j<n; j++){
									if(i==j) null_matrix_column[j]=1.0;
									else null_matrix_column[j]=0.0;
								}
							}
							{
								int info;	/* No need to check info afterwards when job==10000 */
								int job=10000;
								AIR_dqrsl(s2, n, m, qraux, null_matrix_column, job, null_matrix_column, NULL, NULL, NULL, NULL, &info);
							}
							{
								unsigned int j;
								
								for(j=0; j<n; j++){
									s[j][i]=null_matrix_column[j]*sqrt_n;
								}
							}
						}
					}
				}
				{
					/* Solve x'qr_matrix =s' for qr_matrix, which is w'
					 *
					 * (w=s/x)
					 * qr_matrix=x'\s'
					 */
					{
						/* Transpose x into s2
						 *
						 * s2=x'
						 */
						 unsigned int j;
						 
						 for(j=0; j<n; j++){
						 
							unsigned int i;
							
							for(i=0; i<m; i++){
							
								s2[i][j]=x[j][i];
							}
						 }
					}
					{
						/* Perform QR decomposition of s2
						 *
						 * [q, r]=qr(s2,0)
						 */
						 double *qraux=constants[15];		/* vector of length m */
						 
						 AIR_dqrdc(s2, n, m, NULL, NULL, 0, qraux);
						 {
						 	/* Solve for each column of w
						 	 *
						 	 * s[i]'=x'*w[i]
						 	 */
							double *null_matrix_column=constants[16];	/* vector of length n */
							unsigned int i;
							for(i=0; i<m; i++){
								{
									unsigned int j;
									for(j=0; j<n; j++){
										null_matrix_column[j]=s[j][i];
									}
								}
								{
									int info;
									int job=1100;
									AIR_dqrsl(s2, n, m, qraux, null_matrix_column, job, NULL, null_matrix_column, qr_matrix[i], NULL, NULL, &info);
									if(info>-1){
#ifndef USING_R
										printf("singular matrix in ica_MLCf10--use svd to reduce input matrix dimensionality\n");
#else
										REprintf("singular matrix in ica_MLCf10--use svd to reduce input matrix dimensionality\n");
#endif
									}
								}
							}
						}
						
					}
				}
			}
			{
				/* Transpose qr_matrix into w
				 *
				 * qr_matrix=w'
				 */
				unsigned int j;
						
				for(j=0; j<m; j++){
				
					unsigned int i;
					
					for(i=0; i<m; i++){
					
						w[j][i]=qr_matrix[i][j];
					}
				}
			}
		}
		else{
			/* Compute the sources
			 *
			 * s=w*x;
			 */
			AIR_dgemm('n', 'n', m, n, m, 1.0, w, x, 0.0, s);
		}
		/* Compute det(w) (=det(qr_matrix))
		 *
		 * determinant=abs(det(qr_matrix))
		 * f=f-n*log(abs(determinant))
		 */
		{
			unsigned int *ipvt=constants[5]; 	/* length m */
			
			AIR_dgefa(qr_matrix, m, ipvt); 		/* omitting check for singular matrix */
			{
				double determinant=AIR_ddet(qr_matrix, m, ipvt);
				determinant=log(fabs(determinant));
				*logdeterminant=determinant;
				f-=n*determinant;
			}
		}
	}
	/* At this point, s and w are both valid, but gaussian components of both
	 *  would need to be remixed for sources to correspond to PCA
	 *
	 * However, even without remixing, the likelihood and derivatives will be valid
	 */
	{
		f+=subgc*n*0.5*(log_pi+1.0);		/* sub-Gaussian */
		f+=supergc*n*log_2;					/* super-Gaussian */
		f+=(m-ngc)*n*(log_2pi+1.0)/2.0;		/* Gaussian */
		{
			unsigned int j;
			
			for(j=0;j<n;j++){
			
				unsigned int i;
				
				for(i=0; i<subgc; i++){
				
					/* sub-Gaussian
					 */
					f+=s[j][i]*s[j][i];
					f-=logcosh(s[j][i]*sqrt_2);
				}
				for(;i<ngc;i++){
				
					/* super-Gaussian
					 */
					f+=logcosh(s[j][i]*pi/2.0);
				}
				/* Gaussian contribution has already been incorporated
				 *  above as (m-ngc)*n*1.0/2.0 in f+=(m-ngc)*n*(log_2pi+1.0)/2.0;
				 *
				 * f+=s[j][i]*s[j][i]/2.0;
				 */
				
			}
		}
		if(dw_temp==NULL){
		
			/* Remix Gaussian sources by variance to give PCA
			 *  in the process, w is destroyed
			 */
			double *variance=constants[19];
			{
				/* Compute a=inv(w)
				 *
				 * a=inv(w)
				 * w=a
				 */
				unsigned int *ipvt=constants[5]; 	/* length m */

				if(AIR_dgefa(w,m,ipvt)!=m){
					/* This should not happen since singular matrices should already excluded
					 * have been excluded by the calling routine
					 *
					 */
					/* printf("%s: %c: Error in computing inv(W')\n",__FILE__,__LINE__); */
				}
				{
					double *work=constants[7];	/* length m */

					AIR_dgedi(w,m,ipvt,work);
				}
			}
			{
				double **a=w;
				unsigned int k;
				
				for(k=0; k<ngc; k++){
					/* Compute variance associated with kth source row of S
					 *
					 * for k=1:ngc
					 *   variance(k)=0;
					 *   for j=1:n
					 *     for i=1:m
					 *      variance(k)=variance(k)+(a(i,k)*s(k,j))^2
					 *     end
					 *   end
					 * end
					 */
					double variance_k=0;
					unsigned int j;
					
					for(j=0; j<n; j++){
					
						unsigned int i;
						
						for(i=0; i<m; i++){
							/* variance_k=variance_k+(w(i,k)*s(k,j))^2
							 */
							double temp=a[k][i]*s[j][k];
							variance_k+=temp*temp;
						}
					}
					variance[k]=variance_k;
				}
				if(gc==1){
					/* Compute variance of the sole Gaussian source
					 *
					 * k=ngc+1
					 * variance(k)=0;
					 * for j=1:n
					 *  for i=1:m
					 *    variance(k)=variance(k)+(a(i,k)*s(k,j))^2
					 *  end
					 * end
					 */
					  
					double variance_k=0;
					
					unsigned int j;
					
					for(j=0; j<n; j++){
					
						unsigned int i;
						
						for(i=0; i<m; i++){
							/* variance_k=variance_k+(a(i,k)*s(k,j))^2
							 */
							double temp=a[k][i]*s[j][k];
							variance_k+=temp*temp;
						}
					}
					variance[k]=variance_k;
				}
			}
			if(gc>1){
				/* Need to remix Gaussian components using SVD */
				{
					/* Drop non-Gaussian columns of a
					 *
					 * asmall=a(:,ngc+1:end)
					 *
					 */
					double **asmall=w+ngc;
					double **v=asmall;
					{
						/* Compute SVD of asmall, destroying asmall in the process
						 *
						 * [u_null,d,v]=svd(asmall);
						 */
						double *e=constants[15];			/* double vector of length min(M, N) (i.e., of length M) */
						double *d=constants[16];			/* vector of length min(n+1, m) (i.e., of lengh m) */
						double *svd_work=constants[17];		/* vector of length n */

						unsigned int info;
						
						AIR_dsvdc(asmall, m, gc, d, e, NULL, v, svd_work, (unsigned int) 0, (unsigned int) 1, &info);
					
						{
							/* Retrieve variance from d 
							 *
							 * for k=(ngc+1):m
							 *   variance(k)=d(k-ngc, k-ngc)*d(k-ngc, k-ngc);
							 * end
							 */
							unsigned int k;
							
							for(k=ngc; k<m; k++){
							
								variance[k]=d[k-ngc]*d[k-ngc]*n;
							}
						}
					}
					{
						double *column=constants[16];		/* vector of lenght m */
						unsigned int j;
						/* Remix Gaussian sources
						 *
						 * s((ngc+1):end,:)=v'*s((ngc+1):end,:)
						 */
						
						for(j=0; j<n; j++){
						
							/* Copy Gaussian rows of column j of s into column */
							{
								/*
								 * column=s(ngc+1:end,j);
								 */
								unsigned int i;
								for(i=ngc; i<m; i++){
									
									column[i-ngc]=s[j][i];
								}
							}
							{
								/*
								 * s(ngc+1:end, j)=v'*column=v'*s(ngc+1:end, j)
								 */
								unsigned int k;
								
								for(k=0; k<gc; k++){
								
									double temp=0;
									
									unsigned int i;
									
									for(i=0; i<gc; i++){
									
										temp+=v[k][i]*column[i];
									}
									s[j][k+ngc]=temp;
								}
							}
						}
					}
				}
			}
			/* At this point, s is valid and its Gaussian components
			 *  have been remixed to correspond to PCA with 
			 *  sources in descending order by variance
			 *
			 * The values stored in variance are all correct and
			 *  again correspond to PCA for the Gaussian components
			 *
			 * The value in w is not valid
			 */
			return f;
		}
	
		/* Compute gradient */
		{
			double **dw=constants[6];
			/* Compute dW=inv(W') */
			{

				/* dW=W' */
				{
					unsigned int j;
					
					for(j=0;j<m;j++){
					
						unsigned int i;
						
						for(i=0;i<m;i++){
							dw[i][j]=w[j][i];
						}
					}
				}
				/* dW=inv(W') */
				{
					unsigned int *ipvt=constants[5]; 	/* length m */

					if(AIR_dgefa(dw,m,ipvt)!=m){
						/*
						 * This should not happen since singular matrices are already excluded
						 * printf("%s: %c: Error in computing inv(W')\n",__FILE__,__LINE__);
						 */
					}
					{
						double *work=constants[7];	/* length m */

						AIR_dgedi(dw,m,ipvt,work);
					}
				}
			}
			{
				unsigned int j;
				for(j=0;j<n;j++){
				
					unsigned int i;
					
					for(i=0; i<subgc; i++){
					
						/* Sub-Gaussian */
						s[j][i]=2.0*s[j][i]-sqrt(2.0)*tanh(sqrt(2.0)*s[j][i]);				
					}
					for(;i<ngc;i++){
					
						/* Super-Gaussian */
						s[j][i]=0.5*pi*tanh(s[j][i]*0.5*pi);				
					}
					
					/* Nothing to do for Gaussian sources */
				}
			}
			/* Compute dW=F(S)*X'-N*inv(W') */
			AIR_dgemm('n','t',m,m,n,1.0,s,x,-1.0*n,dw);
					
			
			/* Copy results into dw_temp */
			{
				double *temp=dw_temp;
				unsigned int j;
				
				for(j=0;j<m;j++){
				
					unsigned int i;
					
					for(i=0; i<ngc; i++){
					
						*temp++=dw[j][i];
					}
				}
			}
		}
		return f;
	}
}

