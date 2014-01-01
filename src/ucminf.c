/* Copyright 2002-2010 Roger P. Woods, M.D. */
/* Modified: 09/25/10 */
 
#include "AIR.h"

static void dcopy(const unsigned int n, const double *dx, const int incx, double *dy, const int incy){

	if(n==0) return;

	if(incx<0) dx+=(n-1)*(-incx);
	if(incy<0) dy+=(n-1)*(-incy);
	
	{
		unsigned int i;
		for(i=0;i<n;i++, dx+=incx, dy+=incy){
			*dy=*dx;
		}
	}
}
 
static void dspr2(const unsigned int n, const double *x, const double *y, double *ap){

	/*
	 * Limited port of the FORTRAN version with:
	 *	UPLO='L'
	 *	incx=1;
	 *	incy=1;
	 * 	alpha=1.0
	 *
	 *	A:=ALPHA*X*Y' +ALPHA*Y*X'+A
	 *	so A:=X*Y'+Y*X'
	 */
	
	if(n==0) return;
	
	/* Form A when lower triangle is stored in AP */
	{
		unsigned long int kk=0;
		unsigned int j;
		for(j=0;j<n;kk+=n-j,j++){ // order of kk+=n-j and j++ is important
			if((x[j]!=0.0) || (y[j]!=0.0)){
				double temp1=/* alpha* */ y[j];
				double temp2=/* alpha * */ x[j];
				unsigned long int k=kk;
				unsigned int i;
				for(i=j;i<n;i++,k++){
					ap[k]+=x[i]*temp1+y[i]*temp2;
				}
			}
		}
		return;
	}
}
 
static void sline(double (*fdf)(const unsigned int, const double *, double *, void **), const unsigned int n, double *x, const double f, double *g, void **constants, const double *h, double *w, double *alpha, double *fn, double *slopes, unsigned int *neval)
{

	/* Default return values */
	unsigned int meval=*neval;
	*alpha=0.0;
	*fn=f;
	*neval=0;
	
	/* Get initial slope and check descent direction */
	slopes[0]=AIR_ddot(n,g,h);
	slopes[1]=slopes[0];
	if(slopes[0]>=0.0) return;
	
	/* Split work space and finish initialization */
	{
		double fi0=f;
		double sl0=.05*slopes[0];
		double slthr=.995*slopes[0];
		AIR_Boolean ok=FALSE;
		{
			double xfd[3][3];
			xfd[0][0]=0.0;
			xfd[0][1]=f;
			xfd[0][2]=slopes[0];
			{
				double b=1.0;
				for(;;){
					/* Evaluate at x+b*h */
					xfd[1][0]=b;
					dcopy(n,x,1,w,1);
					AIR_daxpy(n,b,h,w);
					xfd[1][1]=fdf(n,w,w+n,constants);
					*neval=*neval+1;
					xfd[1][2]=AIR_ddot(n,w+n,h);
					if(b==1.0) slopes[1]=xfd[1][2];
					if(xfd[1][1]<=fi0+sl0*xfd[1][0]){

						/* New lower bound */
						if(xfd[1][2]<=fabs(slthr)){

							ok=TRUE;
							*alpha=xfd[1][0];
							*fn=xfd[1][1];
							slopes[1]=xfd[1][2];
							dcopy(n,w+n,1,g,1);
							if((b<2.0) && (xfd[1][2]<slthr)){

								/* Expand */
								dcopy(3,xfd[1],1,xfd[0],1);
								b=2.0;
							}
							else break;
						}
						else break;
					}
					else break;
				}
				{
					double d=xfd[1][0]-xfd[0][0];
					for(;;){
						if(ok || (*neval==meval)) return;
						
						/* Refine interval. Min of quadratic interpolator */
						{
							double c=xfd[1][1]-xfd[0][1]-d*xfd[0][2];
							if(c>1e-15*n*xfd[1][0]){
								/* Minimizer in interval */
								double a=xfd[0][0]-.5*xfd[0][2]*(d*d/c);
								d*=0.1;
								
								xfd[2][0]=xfd[0][0]+d;
								if(a>xfd[2][0]) xfd[2][0]=a;
								if(xfd[1][0]-d<xfd[2][0]) xfd[2][0]=xfd[1][0]-d;
							}
							else{
								xfd[2][0]=.5*(xfd[0][0]+xfd[1][0]);
							}
							dcopy(n,x,1,w,1);
							AIR_daxpy(n,xfd[2][0],h,w);
							xfd[2][1]=fdf(n,w,w+n,constants);
							*neval=*neval+1;
							xfd[2][2]=AIR_ddot(n,w+n,h);
							if(xfd[2][1]< fi0+sl0*xfd[2][0]){
								/* New lower bound */
								ok=TRUE;
								*alpha=xfd[2][0];
								*fn=xfd[2][1];
								slopes[1]=xfd[2][2];
								dcopy(n,w+n,1,g,1);
								dcopy(3,xfd[2],1,xfd[0],1);
							}
							else{
								dcopy(3,xfd[2],1,xfd[1],1);
							}
							/* Check convergence */
							d=xfd[1][0]-xfd[0][0];
							ok=ok && (fabs(xfd[2][2])<=fabs(slthr));
							ok=ok || (d<=0.0);
						}
					}
				}
			}
		}
	}
}
 
static void dspmv(const unsigned int n, const double alpha, const double *ap, const double *x, double *y){
	
	/*
	 * Limited port of the Fortran version with:
	 *	uplo='L'
	 *	incx=1
	 * 	incy=1
	 *	beta=0
	 *
	 * 	y:=alpha*A*x+beta*y
	 *	so, y=alpha*A*x
	 */
	 	
	/* Quick return if possible */
	if(n==0) return;
	
	/* Set up the start points in X and Y */
	
	/* First form y:=beta*y */
	{
		unsigned int i;
		for(i=0;i<n;i++){
			y[i]=0.0;
		}
	}
	if(alpha==0.0) return;
	{
		unsigned long int kk=0;

		/* Form y when AP contains the lower triangle */
		
		unsigned int j;
		for(j=0;j<n;kk+=n-j,j++){ // order of kk+=n-j and j++ is important
			double temp1=alpha*x[j];
			double temp2=0.0;
			y[j]+=temp1*ap[kk];
			{
				unsigned long int k=kk+1;
				unsigned int i;
				for(i=j+1;i<n;i++,k++){
					y[i]+=temp1*ap[k];
					temp2+=ap[k]*x[i];
				}
				y[j]+=alpha*temp2;
			}
		}
	}
}

/*
static void dspr(const unsigned int n, const double alpha, const double *x, double *ap){

	if((n==0) || (alpha==0.0)) return;
				

	{
		unsigned long int kk=0;
		unsigned int j;
		for(j=0;j<n;kk+=n-j,j++){ // order of kk+=n-j and j++ is important
			if(x[j]!=0.0){
				double temp=alpha*x[j];
				unsigned long int k=kk;
				unsigned int i;
				for(i=j;i<n;i++,k++){
					ap[k]+=x[i]*temp;
				}
			}
		}
	}
	return;
}
*/

/*
static unsigned int spchol(const unsigned int n, double *a){
	
	unsigned int kk=0;
	unsigned int k;
	for(k=0;k<n;k++){
		if(a[kk]<=0.0){
			return k;
		}
		a[kk]=sqrt(a[kk]);
		if(k<n-1){
			unsigned int nk=n-k;
			AIR_dscal(nk,1.0/a[kk],a+kk+1);
			{
				unsigned int kn=kk+nk+1;
				dspr(nk,-1.0,a+kk+1,a+kn);
				kk=kn;
			}
		}
	}
	return n;
}
*/

 
static AIR_Boolean chkdfn(double (*fdf)(const unsigned int, const double *, double *, void **), const unsigned int n, double *x, void **constants, double stepl, double *diff, unsigned int *indx, double *g, double *g1){
 	/* Initialize */
 	{
 		unsigned int i;
 		for(i=0;i<4;i++){
 			diff[i]=0.0;
 		}
 	}
 	{
 		unsigned int i;
 		for(i=0;i<3;i++){
 			indx[i]=0;
 		}
 	}
 	{
 		double f=fdf(n,x,g,constants);
 		
 		/* Run through components of x */
 		{
 			unsigned int i;
 			
 			for(i=0;i<n;i++){
 				if(fabs(g[i])>diff[0]) diff[0]=fabs(g[i]);
 				{
 					double ab,af;
 					double xi=x[i];
 					
 					/* Forward */
 					x[i]=xi+stepl;
 					{
 						double h=x[i]-xi;
 						if(h==0.0) return FALSE;
						{
							double f1=fdf(n,x,g1,constants);
							af=(f1-f)/h;
							{
								double er=af-g[i];
								if(fabs(er)>fabs(diff[1])){
									diff[1]=er;
									indx[0]=i;
								}
							}
						}
 					}
 					/* Back */
 					x[i]=xi=0.5*stepl;
 					{
 						double h=x[i]-xi;
 						if(h==0.0) return FALSE;
						{
							double f1=fdf(n,x,g1,constants);
							ab=(f1-f)/h;
							{
								double er=ab-g[i];
								if(fabs(er) > fabs(diff[2])){
									diff[2]=er;
									indx[1]=i;
								}
							}
						}
 					}
 					/* Extrapolated */
 					{
						double ae=(2.0*ab+af)/3.0;
						{
							double er=ae-g[i];
							if(fabs(er) > fabs(diff[3])){
								diff[3]=er;
								indx[2]=i;
							}
						}
 					}
 					/* Restore */
 					x[i]=xi;
 				}
 			}
 			return TRUE;
 		}
 	}
}

int AIR_ucminf(double (*fdf)(const unsigned int, const double *, double *, void **), const unsigned int n, double *x, double *dx, void **constants, double *eps, unsigned int *maxfun, double *w, const unsigned long int iw, const int icontr){

	/*
	 *	This is a C port of the Fortran program UCMINF written by Hans Bruun Nielsen
	 *	See: "UCMINF--AN ALGORITHM FOR UNCONSTRAINED, NONLINEAR OPTIMIZATION, Report
	 *	IMM-REP-2000-18, Department of Mathematical Modelling, Technical University of
	 *	Denmark, December 2000
	 *
	 *	The port is limited in that the feature for initializing D0 is omitted
	 *	In addition a void ** is added to carry constants (such as allocated memory) into the
	 *	function being minimized
	 *
	 *	Also modified is the fact that the value ICONTR is not changed--instead the output
	 *	value that it would normally get is the return value of the function
	 *
	 * 	On Entry:
	 *		FDF is a pointer to a function that returns the functions value as a double
	 *		FDF must be of the form F=FDF(N,X,DX)
	 *
	 *		N is the number of constants used by FDF
	 *
	 *		X is an array for the N constants used by FDF
	 *
	 *		DX a pointer to the gradient
	 *
	 *		If ICONTR is >0, the initial values of X and DX are updated iteratively
	 *		If ICONTR is <0, the Jacobian is checked at X and the gradient is checked
	 *		by finite differences with steplength DX and neither X nor DX are changed
	 *
	 *		EPS is a real array with 2 elements. It is used only when ICONTR >0. The
	 *		algorithms stops based on these values. EPS values are unchanged 
	 *
	 *		MAXFUN is an integer than is only used when ICONTR>0. It is the maximum number
	 *		of calls to FDF. This value is modified
	 *
	 *		W is an array with IW elements. Its entry values are NOT used 
	 *		INIMPLEMENTEDE IS THE FEATURE THAT when ICONTR>2 in 
	 *		which case it should hold an approximation to the lower triangle of the inverse
	 *		of H(x) stored columnwise
	 *
	 *		IW must be at least
	 *			max(N*(N+13)/2,7) if ICONTR<=2
	 *			N*max(N+1,(n+13)/2) if ICONTR>2--NOT IMPLEMENTED
	 *
	 *		ICONTR is an integer that controls computation
	 *			ICONTR<=0: Checks gradient only
	 *			ICONTR>0: Start minimization with:
	 *				D0=I
	 *				ORIGINAL VERSION ALLOWS D0 INITIALIZATION IN W--NOT IMPLEMENTED HERE
	 *
	 *		On Exit:
	 *
	 *			Return value gives performance information:
	 *				1: Stopped by small gradient
	 *				2: Stopped by small step
	 *				3: Stopped because too many iterations (MAXFUN)
	 *				4: Stopped by zero step from line search
	 *				-2: Non-start due to N==0
	 *				-4: Non-start due to fabs(dx) too small
	 *				-5: Non-start due to eps<=0
	 *				-6: Non-start due to MAXFUN=0
	 *				-7: Non-start due to given D not positive definite
	 *				-8: Non-start due to IW too small
	 *
	 *			X contains the computed solution if ICONTR>0
	 *
	 *			DX contains the final gradient if ICONTR>0
	 *
	 *			MAXFUN contains the number of calls made to FDF
	 *
	 *			W[0] contains the maximum element in fabs(F'(x))
	 *			W[1],W[4] delta F and j F
	 *			W[2],W[5] delta B and j B
	 *			W[3],W[6] delta E and j E
	 *
	 *			In case of error, W[4],W[5],W[6] indexes erroneous gradient component
	 */
	 
	AIR_Boolean optim=FALSE;
	//AIR_Boolean dgivn=FALSE;
	AIR_Boolean usedel=FALSE;
	
	unsigned long int nn;
			
	/* What job ? */
	if(icontr>0) optim=TRUE;
	//if(icontr>2) dgivn=TRUE;
	
	/* Simple checks */
	nn=(n*(n+1))/2;
	
	if(n==0) return -2;
	if(optim){
		if(*dx<=0.0) return -4;
		if ((eps[0]<=0.0) || eps[1]<=0.0) return -5;
		if (*maxfun==0) return -6;
	}
	else if (*dx==0.0) return -4;
	else{
		unsigned long int t1=(n*(n+13))/2;
		//unsigned int t2=nn*2;
		
		if(7>t1) t1=7;
		//if(nn+5*n>t2) t2=nn+5*n;
		
		if(iw<t1 /* || (dgivn && (iw<t2))*/ ) return -8;
	}
	if(!optim){
		/* Check gradient */
		unsigned int gn=4;
		unsigned int hn=gn+n;
		{
			unsigned int indx[3];
			if (!chkdfn(fdf,n,x,constants,*dx,w,indx,w+gn,w+hn)) return -4; // DX is too small
			{
				unsigned int i;
				for(i=0;i<3;i++){
					w[i+4]=(double)(indx[i]);
				}
			}
		}
		return 0;
	}

	/* Optimize. Split workspace */
	{
		double *gp=w+n;
		double *gn=gp+n;
		double *hn=gn+n;
		double *dn=hn+n;
		double *wn=dn+nn;
		
		/*
		if(dgivn){
			dcopy(nn,dn,1,w,1);
			dcopy(nn,w+4*n,1,w+4*n+2*nn,1);
			{
				unsigned int temp=spchol(n,w+4*n+2*nn);
				if(temp!=n){
					printf("spchol error at %u\n",temp);
					return -7;
				}
			}
			dcopy(nn,w,-1,dn,-1);
			usedel=FALSE;
		}*/
		//else{

			/* Initialize inverse Hessian to unit matrix */
			{
				unsigned long int i;
				for(i=3*n;i<=4*n+nn;i++){
					w[i]=0.0;
				}
			}
			{
				unsigned int ii=4*n;
				unsigned int di=n;
				unsigned int i;
				for(i=0;i<n;i++){
					w[ii]=1.0;
					ii+=di;
					di-=1;
				}
			}
			usedel=TRUE;

		//}
		/* First call of FDF */
		{
			unsigned int neval=1;
			double fx=fdf(n,x,gn,constants);
			{
				double nmh=0.0;
				double nmx=AIR_dnrm2(n,x);
				double nmg=fabs(*(gn+ AIR_idamax(n,gn)));
				if(nmg<=eps[0]){
					*maxfun=neval;
					w[0]=fx;
					w[1]=nmg;
					w[2]=nmh;
					return 1;
				}
				for(;;){
					/* Copy current x and gradient and get new step */
					dcopy(n,x,1,w,1);
					dcopy(n,gn,1,gp,1);
					dspmv(n,-1.0,dn,gn,hn);
					{
						/* Adjust step length to trust region */
						AIR_Boolean redu=FALSE;
						nmh=AIR_dnrm2(n,hn);
	
						if(nmh<=eps[1]*(eps[1]+nmx)){
							*maxfun=neval;
							w[0]=fx;
							w[1]=nmg;
							w[2]=nmh;
							return 2;
						}
						if((nmh>*dx) || usedel){
							redu=TRUE;
							AIR_dscal(n,*dx/nmh,hn);
	
							nmh=*dx;
							usedel=FALSE;
						}
						/* Line search */
						{
							unsigned int meval=5;
							double a;
							double fxn;
							double sl[2];

							sline(fdf,n,x,fx,gn,constants,hn,wn,&a,&fxn,sl,&meval);
							
							if(a==0.0){
								*maxfun=neval;
								w[0]=fx;
								w[1]=nmg;
								w[2]=nmh;
								return 1;
							}
							/* Update neval, x, f(x) and ||g|| */
							neval=neval+meval;
							nmg=fabs(*(gn+AIR_idamax(n,gn)));
							fx=fxn;
	
							AIR_daxpy(n,a,hn,x);

							nmx=AIR_dnrm2(n,x);
							AIR_daxpy(n,-1.0,x,w);
							nmh=AIR_dnrm2(n,w);
							/* Update trust region */
							if(a<1.0){
								/* Reduce Delta */
								*dx*=.35;
							}
							else if(redu && (sl[1]<.7*sl[0])){
								/* Increase Delta */
								*dx*=3.0;
							}
							/* Update of inverse Hessian by BFGS */
							AIR_dscal(n,-1.0,w);
							AIR_daxpy(n,-1.0,gn,gp);
							{
								double yh=-AIR_ddot(n,gp,w);
								if(yh>1e-8*nmh*AIR_dnrm2(n,gp)){
									dspmv(n,-1.0,dn,gp,hn);
									{
										double yv=-AIR_ddot(n,gp,hn);
										a=(1.0+yv/yh)/yh;
									}
									AIR_dscal(n,-1.0/yh,hn);
									AIR_daxpy(n,.5*a,w,hn);
									dspr2(n,w,hn,dn);
								}
								/* Check stopping criteria */
								{
									double thrx=eps[1]*(eps[1]+nmx);
									if(thrx>*dx){
										*dx=thrx;
									}
									{
										int icontr2=0;
										if(neval>=*maxfun) icontr2=3;
										if(nmh<=thrx) icontr2=2;
										if(nmg<=eps[0]) icontr2=1;
										
										if(icontr2!=0){
											*maxfun=neval;
											w[0]=fx;
											w[1]=nmg;
											w[2]=nmh;
											return icontr2;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
