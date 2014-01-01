/* Copyright 1995-2011 Roger P. Woods, M.D. */
/* Modified 3/13/11 */

/*
 * double dnrm2()
 *
 * This routine is based on the SLATEC (BLAS) Fortran routine with
 *  the same name. The Fortran version is in the public domain.
 *
 * Extensive modifications were required to convert this into readable
 * C code
 */


#include "AIR.h"

#define CUTLO 4.441E-16
#define CUTHI 1.304E19

double AIR_dnrm2(const unsigned int n, const double *x)

{
	if(n==0) return 0.0;
    {
        double	sum=0.0;
        double	xmax;
        unsigned int i;
        
        /*Phase I, Sum is zero*/
        for(i=0;i<n;i++){
            if(x[i]!=0.0) break;
        }

        if(i<n){
            /*Prepare for phase II*/
            xmax=fabs(x[i]);
        }
        else{
            /*Finished, return 0.0*/
            return 0.0;
        }

        /*Phase II, SUM is small*/
        for(;i<n;i++){
            if(fabs(x[i])>CUTLO) break;
            if(fabs(x[i])>xmax){
                sum=1.0+sum*(xmax/x[i])*(xmax/x[i]);
                xmax=fabs(x[i]);
            }
            else{
                sum+=(x[i]/xmax)*(x[i]/xmax);
            }
        }
        {
            double hitest;
            
            if(i<n){
                /*Prepare for phase III*/
                sum=(sum*xmax)*xmax;
                hitest=CUTHI/n;
            }
            else{
                /*Finished, return scaled sum*/
                return xmax*sqrt(sum);
            }

            /*Phase III, SUM is mid-range, no scaling*/
            for(;i<n;i++){
                if(fabs(x[i])>=hitest) break;
                sum+=x[i]*x[i];
            }
        }
        if(i<n){
            /*Prepare for phase IV*/
            sum=(sum/x[i])/x[i];
            xmax=fabs(x[i]);
        }
        else{
            /*Finished, return sqrt of unscaled sum*/
            return sqrt(sum);
        }

        /*Phase IV, SUM is large*/
        for(;i<n;i++){
            if(fabs(x[i])>xmax){
                sum=1.0+sum*(xmax/x[i])*(xmax/x[i]);
                xmax=fabs(x[i]);
            }
            else{
                sum+=(x[i]/xmax)*(x[i]/xmax);
            }
        }

        /*Compute square root and adjust for scaling*/
        return xmax*sqrt(sum);
    }
}
