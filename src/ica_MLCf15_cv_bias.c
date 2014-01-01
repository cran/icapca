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

double AIR_ica_MLCf15_cv_bias(/*@unused@*/ const unsigned int unused, const double *w_parameters, double *dw_temp, void **constants, unsigned int new_elements, double **x2){
	
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
	//const unsigned int gc=m-subgc-supergc;				/* Number of Gaussian components */

	const unsigned int ngc=subgc+supergc; /* Number of non-gaussian components */
	
	double *logdeterminant=(double *)constants[11];		/* log of absolute value of determinant is stored here */
	
	double f=0.0;
		
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
		{
			/* Compute the sources
			 *
			 * s=w*x;
			 */
			AIR_dgemm('n', 'n', m, new_elements, m, 1.0, w, x2, 0.0, s);
		}
		/* Compute det(w) (=det(qr_matrix))
		 *
		 * determinant=abs(det(qr_matrix))
		 * f=f-new_elements*log(abs(determinant))
		 */
		{
			unsigned int *ipvt=constants[5]; 	/* length m */
			
			AIR_dgefa(qr_matrix, m, ipvt); 		/* omitting check for singular matrix */
			{
				double determinant=AIR_ddet(qr_matrix, m, ipvt);
				determinant=log(fabs(determinant));
				*logdeterminant=determinant;
				f-=new_elements*determinant;
			}
		}
	}
	/* At this point, s and w are both valid, but gaussian components of both
	 *  would need to be remixed for sources to correspond to PCA
	 *
	 * However, even without remixing, the likelihood and derivatives will be valid
	 */
	{
		f+=subgc*new_elements*0.5*(log_pi+1.0);		/* sub-Gaussian */
		f+=supergc*new_elements*log_2;					/* super-Gaussian */
		f+=(m-ngc)*new_elements*(log_2pi)/2.0;		/* Gaussian */
		{
			unsigned int j;
				
			for(j=0;j<new_elements;j++){
				
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
	
				for(; i<m; i++){
					f+=s[j][i]*s[j][i]/2.0;
				}
					
			}
		}
		return f;
	}
}

