/* Copyright 1996-2011 Roger P. Woods, M.D. */
/* Modified 4/24/2011 */

/*
 * singular value decomposition of an n by p rectangular matrix
 *
 * Note that u is not referenced if joba==0.
 * If n<=p or joba==2, u may be identified with x in the subroutine call.
 *
 * Note that v is not referenced if jobb==0.
 * If p<=n, v may be identified with x in the subroutine call.
 *
 * The input matrix x[p][n] is destroyed.
 * The first min(n,p) entries in s are singular values
 * e is usually zeros
 * u has the singular vectors x[n][n] or x[p][n] depending on joba
 * v has the other singular vectors x[p][p]
 * work is scratch space
 * 
 * If joba==0, no singular vectors are computed
 * If joba==1, the n left singular vectors are returned in u
 * If joba>=2 the first min(n,p) singular vectors are returned in u
 *
 * If jobb==0, no right singular vectors are returned
 * If jobb==1, the right singular vectors are returned in v
 *
 * On return, *info==0 implies all values are correct
 *
 * X is the N x P matrix (x[p][n]) whose singular value is to be computed
 *
 * S will return the MIN(N,P) (s[n+1] or s[p]) singular values in descending magnitude
 * Note that the array must be N+1 if P>=N+1, even though only N values are returned
 *
 * E (e[p]) will contain zeros usually
 *
 * U:
 *  If joba==0: Unused
 * 	If joba==1: (N X N) (u[n][n])
 * 	If joba>=2: (N X MIN(N,P)) (u[n][n] or u[p][n])
 *  Note: if N<=P or joba==2, X and U can be the same
 *
 * V:
 *	If jobb==0: Unused
 *  Otherwise (P X P) (v[p][p])
 *	Note: if P<=N, X and V can be the same	
 *
 * WORK:
 * 	Needs to be length N (work[n])
 */

#include "AIR.h"
#define MAXIT 60	/* Maximum number of iterations */

AIR_Error AIR_dsvdc(double **x, const unsigned int n, const unsigned int p, double *s, double *e, double **u, double **v, double *work, const unsigned int joba, const AIR_Boolean jobb, unsigned int *info)

{
	unsigned int m;
	
	/*Determine what is to be computed*/
	AIR_Boolean wantu=FALSE;
	AIR_Boolean wantv=FALSE;
	unsigned int ncu=n;
	
	if(joba>1){
		if(p<ncu) ncu=p;
	}
	if(joba!=0) wantu=TRUE;
	if(jobb) wantv=TRUE;

	/*Reduce x to bidiagonal form, storing the diagonal elements in */
	/* s and the super-diagonal elements in e */

	if(info!=NULL) *info=0;
	if(n==0) return 0;
	{
		unsigned int nct=n-1;
		if(p<nct) nct=p;
		{
			unsigned int nrt=0;
			if(p>1){
				nrt=p-2;
				if(n<nrt) nrt=n;
			}
			{
				unsigned int lu=nct;
				if(nrt>lu) lu=nrt;
				if(lu>=1){
				
					unsigned int l;
					
					for(l=0;l<lu;l++){
						if(l+1<=nct){
							/*Compute the transformation for the l-th */
							/*column and place the lth diagonal in s(l) */
							s[l]=AIR_dnrm2(n-l,&x[l][l]);
							if(s[l]!=0.0){
								if(x[l][l]!=0.0){
									if(x[l][l]>0) s[l]=fabs(s[l]);
									else s[l]=-fabs(s[l]);
								}
								AIR_dscal(n-l,1.0/s[l],&x[l][l]);
								x[l][l]+=1.0;
							} /*10*/
							s[l]=-s[l];
						} /*20*/
						if(p>l+1){
						
							unsigned int j;
							
							for(j=l+1;j<p;j++){
								if(l+1<=nct){
									if(s[l]!=0.0){
										/*apply transformation*/
										double t=-AIR_ddot(n-l,&x[l][l],&x[j][l])/x[l][l];
										AIR_daxpy(n-l,t,&x[l][l],&x[j][l]);
									}
								} /*30*/
								/*place lth row of x into e for the subsequent */
								/*calculation of the row transformation		*/
								e[j]=x[j][l];
							} /*40*/
						} /*50*/
						if(wantu && l+1<=nct){
							/*place the transformation in u for subsequent back multiplication*/
							unsigned int i;
							
							for(i=l;i<n;i++){
								u[l][i]=x[l][i];
							}
						}/*70*/
						if(l+1<=nrt){
							/*compute the lth row transformation and place the lth superdiagonal*/
							/*in e[l] 								*/
							e[l]=AIR_dnrm2(p-l-1,&e[l+1]);
							if(e[l]!=0.0){
								if(e[l+1]!=0.0){
									if(e[l+1]>0.0) e[l]=fabs(e[l]);
									else e[l]=-fabs(e[l]);
								}
								AIR_dscal(p-l-1,1.0/e[l],&e[l+1]);
								e[l+1]+=1.0;
							}/*80*/
							e[l]=-e[l];
							if(l+2<=n && e[l]!=0.0){
								/*apply transformation*/
								{
									unsigned int i;
									
									for(i=l+1;i<n;i++){
										work[i]=0.0;
									}
								}
								{
									unsigned int j;
									
									for(j=l+1;j<p;j++){	
										AIR_daxpy(n-l-1,e[j],&x[j][l+1],&work[l+1]);
									}
								}
								{
									unsigned int j;
									
									for(j=l+1;j<p;j++){
										AIR_daxpy(n-l-1,-e[j]/e[l+1],&work[l+1],&x[j][l+1]);
									}
								}
							}/*120*/
							if(wantv){
								/*place transformation in v for subsequent back multiplication*/
								
								unsigned int i;
								
								for(i=l+1;i<p;i++){
									v[l][i]=e[i];
								}
							}/*140*/
						}/*150*/
					}
				} /* 160 */ /*170*/
			}

			
			/*Set up the final bidiagonal matrix or order m*/
			{
				m=n+1;
				if(p<m) m=p;
				if(nct<p) s[nct]=x[nct][nct];
				if(n<m) s[m-1]=0.0;
				if(nrt+1<m) e[nrt]=x[m-1][nrt];
				e[m-1]=0.0;

				/*if required, generate u */
				if(wantu){
					if(ncu>=nct+1){
						
						unsigned int j;
						
						for(j=nct;j<ncu;j++){
							
							unsigned int i;
							
							for(i=0;i<n;i++){
								u[j][i]=0.0;
							}
							u[j][j]=1.0;
						}
					}/*200*/
					if(nct>=1){
					
						unsigned int ll;
						
						for(ll=0;ll<nct;ll++){
						
							unsigned int l=nct-ll-1;
							if(s[l]!=0.0){
								if(ncu>=l+2){
								
									unsigned int j;
									
									for(j=l+1;j<ncu;j++){
										double t=-AIR_ddot(n-l,&u[l][l],&u[j][l])/u[l][l];
										AIR_daxpy(n-l,t,&u[l][l],&u[j][l]);
									}
								}/*220*/
								AIR_dscal(n-l,-1.0,&u[l][l]);
								u[l][l]+=1.0;
								if(l>=1){
								
									unsigned int i;
									
									for(i=0;i<l;i++){
										u[l][i]=0.0;
									}/*230*/
								}/*240*/
							}/*250*/
							else{
								{
									unsigned int i;
									
									for(i=0;i<n;i++){
										u[l][i]=0.0;
									}
								}
								u[l][l]=1.0;
							}/*270*/
						}/*280*/
					}/*290*/
				}/*300*/

				/*if required, generate v*/

				if(wantv){
					
					unsigned int ll;
					
					for(ll=0;ll<p;ll++){
					
						unsigned int l=p-ll-1;
						if(l+1<=nrt){
							if(e[l]!=0.0){
							
								unsigned int j;
								
								for(j=l+1;j<p;j++){
								
									double t=-AIR_ddot(p-l-1,&v[l][l+1],&v[j][l+1])/v[l][l+1];
									AIR_daxpy(p-l-1,t,&v[l][l+1],&v[j][l+1]);
								}
							}
						}/*320*/
						{
							unsigned int i;
							
							for(i=0;i<p;i++){
								v[l][i]=0.0;
							}
						}
						v[l][l]=1.0;
					}/*340*/
				}/*350*/
			}
		}
	}

	/*Main iteration loop for the singular values*/
	{
		unsigned int mm=m;
		unsigned int iter=0;
		while(m!=0){	/*360*/
			if(iter>=MAXIT){
				if(info!=NULL) *info=m;
				return AIR_SVD_FAILURE_ERROR;
			}
			{
				unsigned int kase;
				unsigned int l;
				
				for(l=m-1;l!=0;l--){
					double test=fabs(s[l-1])+fabs(s[l]);
					double ztest=test+fabs(e[l-1]);
					
					if(ztest==test){
						e[l-1]=0.0;
						break;
					}
				}
				if(l==m-1) kase=4; /* Convergence: e[m-2] is negligible */
				else{
					/* At this point, m>1, l<m-1 */
					
					unsigned int ls;
					
					for(ls=m;ls>l;ls--){
					
						double test=0.0;
						if(ls!=m) test+=fabs(e[ls-1]);
						if(ls!=l+1) test+=fabs(e[ls-2]);
						{
							double ztest=test+fabs(s[ls-1]);
							if(ztest==test){
								s[ls-1]=0.0;
								break;
							}
						}
					}
					if(ls==l) kase=3;
					else{
						if(ls==m) kase=1;
						else{
							kase=2;
							l=ls;
						}
					}
				}

				/*Perform the task indicated by kase */
				if(kase==1){
					/*Deflate negligible s[m-1] */
					double f=e[m-2];
					e[m-2]=0.0;
					{
						unsigned int kk;
						
						for(kk=l;kk<m-1;kk++){
						
							unsigned int k=m-2-kk+l;
							double t1=s[k];
							{
								double sn,cs;
								
								AIR_drotg(&t1,&f,&cs,&sn);
								s[k]=t1;
								if(k!=l){
									f=-sn*e[k-1];
									e[k-1]*=cs;
								}
								if(wantv) AIR_drot(p,v[k],v[m-1],cs,sn);
							}
						}
					}
				}
				else if(kase==2){
					/*Split at negligible s[l]*/
					double f=e[l-1];
					e[l-1]=0.0;
					{
						unsigned int k;
						
						for(k=l;k<m;k++){
							double t1=s[k];
							{
								double sn,cs;
								
								AIR_drotg(&t1,&f,&cs,&sn);
								s[k]=t1;
								f=-sn*e[k];
								e[k]*=cs;
								if(wantu) AIR_drot(n,u[k],u[l-1],cs,sn);
							}
						}
					}
				}
				else if(kase==3){
					/*Perform one qr step*/
					/*calculate the shift*/
					double sm;
					double emm1;
					double smm1;
					double sl;
					double el;
					double b,c,g;
					double f;
					{
						double scale=fabs(s[m-1]);
						if(fabs(s[m-2])>scale) scale=fabs(s[m-2]);
						if(fabs(e[m-2])>scale) scale=fabs(e[m-2]);
						if(fabs(s[l])>scale) scale=fabs(s[l]);
						if(fabs(e[l])>scale) scale=fabs(e[l]);
						sm=s[m-1]/scale;
						smm1=s[m-2]/scale;
						emm1=e[m-2]/scale;
						sl=s[l]/scale;
						el=e[l]/scale;
					}
					b=((smm1+sm)*(smm1-sm)+emm1*emm1)/2.0;
					c=(sm*emm1)*(sm*emm1);
					{
						double shift=0.0;
						if(b!=0.0 || c!=0.0){
							shift=sqrt(b*b+c);
							if(b<0.0) shift=-shift;
							shift=c/(b+shift);
						}
						f=(sl+sm)*(sl-sm)-shift;
					}
					g=sl*el;
					
					/*Chase zeros*/
					{
						unsigned int k;
						
						for(k=l;k<m-1;k++){
							{
								double sn,cs;
								
								AIR_drotg(&f,&g,&cs,&sn);
								if(k!=l) e[k-1]=f;
								f=cs*s[k]+sn*e[k];
								e[k]=cs*e[k]-sn*s[k];
								g=sn*s[k+1];
								s[k+1]*=cs;
								if(wantv) AIR_drot(p,v[k],v[k+1],cs,sn);
							}
							{
								double sn,cs;
								
								AIR_drotg(&f,&g,&cs,&sn);
								s[k]=f;
								f=cs*e[k]+sn*s[k+1];
								s[k+1]=-sn*e[k]+cs*s[k+1];
								g=sn*e[k+1];
								e[k+1]*=cs;
								if(wantu && (k+1)<n) AIR_drot(n,u[k],u[k+1],cs,sn);
							}
						}/*560*/
					}
					e[m-2]=f;
					iter++;
				}
				else{
					/*Make the singular value positive*/
					if(s[l]<0.0){
						s[l]=-s[l];
						if(wantv) AIR_dscal(p,-1.0,v[l]);
					}
					/*Order the singular value*/
					while((l+1)!=mm){
						if(s[l]>=s[l+1]) break;
						{
							double t=s[l];
							s[l]=s[l+1];
							s[l+1]=t;
						}
						if(wantv && (l+1)<p) AIR_dswap(p,v[l],v[l+1]);
						if(wantu && (l+1)<n) AIR_dswap(n,u[l],u[l+1]);
						l++;
					}/*600*/
					iter=0;
					m--;
				}
			}
		}
	}/*620*/
	return 0;
}
