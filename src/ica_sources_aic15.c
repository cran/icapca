/* Copyright 2005-2014 Roger P. Woods, M.D. */
/* Modified: 10/19/14 */

/*
 * Main algorithm for ICA minimization with sub-Gaussian, super-Gaussian and Gaussian sources
 *
 * On INPUT
 *
 *		SEED			if >0, initializes a pseudorandom orthonormal rotation of input matrix before solving sources
 *
 *		OFFSET_RANDOM	the number of pseudorandom orthonormal rotations to skip before selecting one
 *
 *		SAMPLE_ORDER[N]	an array of length N containing the numbers 1 to N in a suitably randomized order, only used when inf_crit==3
 *							a non-null sample_order will be used if samples==0 (treated as if samples==n)
 *
 *		SAMPLE_OFFSET	the number of elements in sample_order to skip before sampling, only used when inf_crit==3
 *							ignored if samples==0
 *
 *		SAMPLES			the number of elements to sample from sample_order, only used when inf_crit==3, 
 *							if 0 (or n) every observation is left out in turn
 *							otherwise, only the specified samples from sample_order are left out
 *
 *		CONSTANTS_PTR	a void ***
 *
 *		INF_CRIT		the information criterion to use:
 *							0: No correction
 *							1: Akaike Information Criterion uncorrected (AIC)
 *							2: Bayesian Information Criterion of Schwarz (BIC)
 *							3: N-fold leave one out cross validation
 *							>3: Second order AIC correction (requires no non-Gaussian sources
 *
 *		FOLD			when inf_crit==3, specifies the fold of the cross-validation:
 *							0: uses leave-one-out cross validation, same as fold==n, but computed with more efficiently
 *							1: invalid value
 *							1<fold<=n: fold of cross-validation if (n/fold)*fold==n, otherwise invalid value
 *							n: same as fold==0 but computed using generic implementation
 *
 *		XVAL_EPSILON	value used when inf_crit==3 to screen uncorrected loglikelihoods
 *						for redundant initializations
 *							>0.0	assume loglikelihoods differing by less than this amount
 *									 are redundant
 *							0.0		compute cross validation for all initializations
 *							<0.0	assume initializations producing same loglikelihood when
 *									 formatted as sprintf '%e' are redundant
 *
 *		RANGES			should be a pointer of an unsigned int array
 *							specifying the allowable numbers of each component
 *							RANGES[0], RANGES[1]	sub-Gaussian minimum and maximum
 *							RANGES[2], RANGES[3]	super-Gaussian minimum and maximum
 *							RANGES[4], RANGES[5]	Gaussian minimum and maximum
 *
 *
 *		EM[N][POINTS]	the POINTS by N matrix to decompose
 *
 *		POINTS			the number of rows in EM
 *
 *		N				the number of columns in EM
 *
 *		HINTED_SUBGAUSSIAN_SOURCES
 *						the number of hinted subgaussian sources in SOURCE_MATRIX
 *
 *		HINTED_SUPERGAUSSIAN_SOURCES
 *						the number of hinted supergaussian sources in SOURCE_MATRIX
 *
 *		HINTED_UNSPECIFIED_SOURCES
 *						the number of hinted unspecified sources in SOURCE_MATRIX
 *
 *		SOURCE_MATRIX[HINTED_SUBGAUSSIAN_SOURCES+HINTED_SUPERGAUSSIAN_SOURCES+HINTED_UNSPECIFIED_SOURCES][N]		
 *						the N by HINTED_SUBGAUSSIAN_SOURCES+HINTED_SUPERGAUSSIAN_SOURCES+HINTED_UNSPECIFIED_SOURCES
 *						matrix of sources
 *						sources must be in the following order:
 *							hinted_subgaussian_sources
 *							hinted_supergaussian_sources
 *							hinted_unspecified_sources
 *						SOURCE_MATRIx can be NULL if HINTED_SUBGAUSSIAN_SOURCES+HINTED_SUPERGAUSSIAN_SOURCES
 *						+HINTED_UNSPECIFIED_SOURCES==0)
 *
 *		*COMPONENTS		the maximum number of components to retain
 *
 *		ICA_S_PTR		a double *** where results will be placed
 *
 *		DESIRED_INITIALIZATION	an unsigned integer
 *							0 for normal computation
 *							the 1 offset index of the desired initialization among permutations
 *
 *		DISTRIBUTION_PTR
 *						an int ** where results will be placed
 *
 *		VARIANCE_PTR	a double ** where results will be placed
 *
 *		PROBABILITY_PTR	a double ** where results will be placed
 *
 *		LOGLIKELIHOOD	a double * where a result will be placed
 *
 *		OUTPREFIX		the basename for postscript graph files
 *						 if NULL, no graphs will be created
 *
 *		DATA_SAVER		a routine to save intermediate ICA results for each possible number
 *						of Gaussian components (see ica_plot.c for the format required).
 *
 *		PLOT_CONSTANTS	array of void pointers that can carry information through this routine
 *						into data_saver without modification
 *
 *		MODEL_OF_INTEREST
 *						if not NULL, must be an int array of length 3, with entries
 *						specifying the number of sources of each category as follows:
 *							MODEL_OF_INTEREST[0] => sub-Gaussian
 *							MODEL_OF_INTEREST[1] => super-Gaussian
 *							MODEL_OF_INTEREST[2] => Gaussian (orthogonal to all other sources)
 *
 *						A value of -1 serves as a wild card
 *						If MODEL_OF_INTEREST is not NULL, the probability of the model will
 *						be returned in *MOI_PROBABILITY
 *
 *		MOI_PROBABILITY	a double * where a result will be placed
 *
 *		LEAVE_ROWS_UNCENTERED
 *						a boolean indicating whether the rows of em should be left uncentered (default is to center)
 *
 *		VERBOSE			a boolean determining whether to print information about each type of model (default is false)
 *
 *		OW				whether overwrite permission is granted for data_saver output (e.g., postscript graph files)
 *
 * On RETURN
 *
 *		CONSTANTS_PTR	call AIR_free_ica_sources_aic15 with this value to free memory allocated by this routine
 *						this includes the memory used to store *ICA_S_PTR, *DISTRIBUTION_PTR, *VARIANCE_PTR and
 *						*PROBABILITY_PTR
 *
 *		EM				the original contents of this matrix will have been destroyed
 *
 *		*ICA_S_PTR[N][*COMPONENTS] 				
 *						the *COMPONENTS by N source matrix S
 *						do not free the associated storage directly--free constants_ptr instead
 *
 *		*DISTRIBUTION_PTR[*COMPONENTS]			
 *						classification of the sources
 *							0 => sub-Gaussian
 *							1 => super-Gaussian
 *							2 => Gaussian (orthogonal to all other sources)
 *						do not free the associated storage directly--free constants_ptr instead
 *
 *		*VARIANCE_PTR[*COMPONENTS]
 *						the proportion of total variance attributable to each source
 *						do not free the associated storage directly--free constants_ptr instead
 *
 *		*PROBABILITY_PTR[*COMPONENTS+1]
 *						relative probabilities of models with different numbers of Gaussian sources
 *						ordered by the number of Gaussian sources, starting with zero
 *						do not free the associated storage directly--free constants_ptr instead
 *
 *		*LOGLIKELIHOOD		the log likelihood of the most likely model
 *
 *		*COMPONENTS			the number of components retained after initial SVD
 *
 *		SOURCE_MATRIX		the input source matrix will have been multiplied by a scalar value (the largest absolute value in the input matrix)
 *
 *		*MOI_PROBABILITY	the probability of the model of interest, if MODEL_OF_INTEREST was not
 *							NULL. Otherwise, this address is not referenced
 *
 *	NOTE:
 *		1. If sources are provided, each source (i.e., each row of SOURCE_MATRIX) should sum to zero
 *		2. Sources need not be orthogonal, but should not be degenerate
 *		3. The mixing matrix, ica_a is freed and not returned because this will not be the expected
 *			matrix due to the SVD done at the outset (U is not computed, so the expected A cannot be
 *			recovered without additional computations
 *		4. If OUTPREFIX is not NULL, calling routines may want to unblock SIGPIPE signals after 
 *			calling this subroutine as follows:
 *				#include <signal.h>
 *				signal(SIGPIPE, SIG_DEF);
 *		5. If all sources are of the same type (e.g., all supergaussian), this routine will
 *		     perform ICA optimization of the best fitting (i.e., only) model only once. Otherwise
 *		     it will be computed more than once to avoid having to store the result when comparing
 *		     different models
 */

#include "AIR.h"
#include <float.h>
#define CONSTANT_COUNT 100
#define NULL_INDEX 10

#define MAXFUN 25000
 
struct hoods{
	double value;
	double logdeterminant;
	unsigned int rank;
	unsigned int subgc;
	unsigned int supergc;
	unsigned int gc;
};

static void ica_pca_free(void *ptr){
#ifndef USING_R
#ifndef USING_MATLAB
	free(ptr);
#else
	mxFree(ptr);
#endif
#endif
}
static void *ica_pca_malloc(size_t n, int size){
#ifndef USING_R
#ifndef USING_MATLAB
    size_t total=n*size;
    if (total==0) return NULL;
	return(malloc(n*size));
#else
	return(mxMalloc(n*size));
#endif
#else
	return(R_alloc(n, size));
#endif
}

static void ica_pca_free_1(void *pixels)
{
	ica_pca_free(pixels);
}
static void ica_pca_free_2(void **pixels)
{
	if(*pixels) ica_pca_free_1(*pixels);
	ica_pca_free(pixels);
}

static double *ica_pca_matrix1(const size_t a)
{
	double *high=(double *)ica_pca_malloc(a, sizeof(double));
	if(!high) return NULL;
	return high;
}
static double **ica_pca_matrix2(const size_t b, const size_t a)
{
	double **high=(double **)ica_pca_malloc(a, sizeof(double *));

	if(!high) return NULL;

	if(b!=0){
		double *low=ica_pca_matrix1(b*a);

		if(!low){
			ica_pca_free(high);
			return NULL;
		}
		{
			double **ptr;
			
			for(ptr=high; ptr<high+a; ptr++, low+=b){
			
				*ptr=low;
			}
		}
	}
	return high;
}
static int compare_hoods(const void *a, const void *b){

	struct hoods hood1=*(struct hoods *)a;
	struct hoods hood2=*(struct hoods *)b;
	
	double diff=hood1.value-hood2.value;
	
	if(diff==0) return 0;
	if(diff<0.0) return 1;
	return -1;
}
static int hoods_by_index(const void *a, const void *b){
	struct hoods hood1=*(struct hoods *)a;
	struct hoods hood2=*(struct hoods *)b;

	if(hood1.rank==hood2.rank) return 0;
	if(hood1.rank>hood2.rank) return 1;
	return -1;
}

static void set_kpvt(const unsigned int dimensions, unsigned int *distribution, unsigned int *kpvt){
	
	unsigned int index=0;
	{
		unsigned int i;
		
		for(i=0; i<dimensions; i++){
			if(distribution[i]==3){
				kpvt[index++]=i;
			}
		}
		for(i=0; i<dimensions; i++){
			if(distribution[i]==0){
				kpvt[index++]=i;
			}
		}
		for(i=0; i<dimensions; i++){
			if(distribution[i]==1){
				kpvt[index++]=i;
			}
		}
		for(i=0; i<dimensions; i++){
			if(distribution[i]==2){
				kpvt[index++]=i;
			}
		}
	}
}

static void count_gaussians_and_nongaussians(const unsigned int dimensions, unsigned int *distribution, unsigned int *subgc, unsigned int *supergc, unsigned int *gc){
	
	unsigned int mysubs=0;
	unsigned int mysupers=0;
	unsigned int mygaussians=0;
	
	unsigned int i;
	
	for(i=0; i<dimensions; i++){
	
		unsigned int temp=distribution[i];
		if(temp==0) mysubs++;
		else if(temp==1) mysupers++;
		else if(temp==2) mygaussians++;
	}
	*subgc=mysubs;
	*supergc=mysupers;
	*gc=mygaussians;
}

static unsigned int increment_index(const unsigned int dimensions, unsigned int *distribution, unsigned int *ranges, unsigned int *zeros, unsigned int *ones, unsigned int *twos){
	/* Returns 1 if successful
	 *  returns 0 when all indices have been considered
	 */ 
	{
		int i;
		for(i=dimensions-1; i>=0; i--){
		
			/* Remove this item from the count */
			switch(distribution[i]){
				case 0:
					(*zeros)--;
					if(*zeros>ranges[1]) continue;
					break;
				case 1:
					(*ones)--;
					if(*ones>ranges[3]) continue;
					break;
				case 2:
					(*twos)--;
					if(*twos>ranges[5]) continue;
					break;
				case NULL_INDEX:
					if(i==0) break;
					continue;
				default:
					break;
			}
			{
				/* Try to increment at this position */
				int success=0;
				if(distribution[i]==NULL_INDEX){
					if(*zeros<ranges[1]){
						unsigned int mandatory_slots=0;
						if(ranges[0]>*zeros+1) mandatory_slots+=ranges[0]-*zeros-1;
						if(ranges[2]>*ones) mandatory_slots+=ranges[2]-*ones;
						if(ranges[4]>*twos) mandatory_slots+=ranges[4]-*twos;
						if(i+mandatory_slots<dimensions){
							distribution[i]=0;
							(*zeros)++;
							success=1;
						}
					}
					if(!success && *ones<ranges[3]){
						unsigned int mandatory_slots=0;
						if(ranges[0]>*zeros) mandatory_slots+=ranges[0]-*zeros;
						if(ranges[2]>*ones+1) mandatory_slots+=ranges[2]-*ones-1;
						if(ranges[4]>*twos) mandatory_slots+=ranges[4]-*twos;
						if(i+mandatory_slots<dimensions){
							distribution[i]=1;
							(*ones)++;
							success=1;
						}
					}
					if(!success && *twos<ranges[5]){
						unsigned int mandatory_slots=0;
						if(ranges[0]>*zeros) mandatory_slots+=ranges[0]-*zeros;
						if(ranges[2]>*ones) mandatory_slots+=ranges[2]-*ones;
						if(ranges[4]>*twos+1) mandatory_slots+=ranges[4]-*twos-1;
						if(i+mandatory_slots<dimensions){
							distribution[i]=2;
							(*twos)++;
							success=1;
						}
					}
				}
				if(!success && distribution[i]==0){
					if(*ones<ranges[3]){
						unsigned int mandatory_slots=0;
						if(ranges[0]>*zeros) mandatory_slots+=ranges[0]-*zeros;
						if(ranges[2]>*ones+1) mandatory_slots+=ranges[2]-*ones-1;
						if(ranges[4]>*twos) mandatory_slots+=ranges[4]-*twos;
						if(i+mandatory_slots<dimensions){
							distribution[i]=1;
							(*ones)++;
							success=1;
						}
					}
					if(!success && *twos<ranges[5]){
						unsigned int mandatory_slots=0;
						if(ranges[0]>*zeros) mandatory_slots+=ranges[0]-*zeros;
						if(ranges[2]>*ones) mandatory_slots+=ranges[2]-*ones;
						if(ranges[4]>*twos+1) mandatory_slots+=ranges[4]-*twos-1;
						if(i+mandatory_slots<dimensions){
							distribution[i]=2;
							(*twos)++;
							success=1;
						}
					}
				}
				if(!success && distribution[i]==1){
					if(*twos<ranges[5]){
						unsigned int mandatory_slots=0;
						if(ranges[0]>*zeros) mandatory_slots+=ranges[0]-*zeros;
						if(ranges[2]>*ones) mandatory_slots+=ranges[2]-*ones;
						if(ranges[4]>*twos+1) mandatory_slots+=ranges[4]-*twos-1;
						if(i+mandatory_slots<dimensions){
							distribution[i]=2;
							(*twos)++;
							success=1;
						}
					}
				}
				if(success){
					unsigned int remaining_dimensions=dimensions;
					while(ranges[4]>*twos){
						distribution[remaining_dimensions-1]=2;
						remaining_dimensions--;
						(*twos)++;
					}
					while(ranges[2]>*ones){
						distribution[remaining_dimensions-1]=1;
						remaining_dimensions--;
						(*ones)++;
					}
					{
						int j;
						for(j=i+1; j<remaining_dimensions; j++){
							/* Insert the smallest available number */
							if(*zeros<ranges[1]){
								distribution[j]=0;
								(*zeros)++;
							}
							else if(*ones<ranges[3]){
								distribution[j]=1;
								(*ones)++;
							}
							else if(*twos<ranges[5]){
								distribution[j]=2;
								(*twos)++;
							}
							else{
#ifndef USING_R
#ifndef USING_MATLAB
								printf("%s: %d: ",__FILE__,__LINE__);
								printf("Coding error in ICA routine\n");
								exit(1);
#endif
#endif
							}
						}
					}
					return(1);
				}
			}
		}
	}
	{
		unsigned int i;
		for(i=0; i<dimensions; i++){
			distribution[i]=NULL_INDEX;
		}
		return(0);
	}
}
static unsigned int count_valid_indices_quickly(const unsigned int dimensions, unsigned int *distribution, unsigned int *ranges){
	// compute dimesions!/(subgc!*supergc!*gc!)
	unsigned int count=0;
	double count2=0.0;
	
	if(ranges[0]+ranges[2]+ranges[4]==dimensions){
	
		unsigned int r0=ranges[0];
		unsigned int r1=ranges[2];
		unsigned int r2=ranges[4];
		
		// Sort so r2>=r1>=r0
		if(r0>r1){
			unsigned int temp=r1;
			r1=r0;
			r0=temp;
		}
		if(r0>r2){
			unsigned int temp=r2;
			r2=r0;
			r0=temp;
		}
		if(r1>r2){
			unsigned int temp=r2;
			r2=r1;
			r1=temp;
		}
		
		count=1;
		/* Protect against overflow */
		count2=1.0;
		{
			unsigned int i;
			for(i=r2+1; i<=dimensions; i++){
				count*=i;
				count2*=i;
			}
		}
		if(count2>(double)UINT_MAX){
			count=UINT_MAX;
		}
		else{
			{
				unsigned int i;
				for(i=r1; i>1; i--){
					count/=i;
					count2/=i;
				}
			}
			{
				unsigned int i;
				for(i=r0; i>1; i--){
					count/=i;
					count2/=i;
				}
			}
		}
	}
	{
		int i;
		for(i=0; i<dimensions; i++){
			distribution[i]=NULL_INDEX;
		}
	}
	return count;
}

static unsigned int count_valid_indices(const unsigned int dimensions, unsigned int *distribution, unsigned int *ranges){
	if(ranges[0]==ranges[1] && ranges[2]==ranges[3] && ranges[4]==ranges[5]){
		return count_valid_indices_quickly(dimensions, distribution, ranges);
	}
	else{
		unsigned int count=0;
		
		unsigned int zeros=0;
		unsigned int ones=0;
		unsigned int twos=0;
		
		unsigned int max_subgc=0;
		unsigned int min_subgc=(unsigned int)-1;
		unsigned int max_supergc=0;
		unsigned int min_supergc=(unsigned int)-1;
		unsigned int max_gc=0;
		unsigned int min_gc=(unsigned int)-1;
						
		{
			int i;
			for(i=0; i<dimensions; i++){
				distribution[i]=NULL_INDEX;
			}
		}
		while(increment_index(dimensions, distribution, ranges, &zeros, &ones, &twos)!=0){
			if(zeros>max_subgc) max_subgc=zeros;
			if(zeros<min_subgc) min_subgc=zeros;
			if(ones>max_supergc) max_supergc=ones;
			if(ones<min_supergc) min_supergc=ones;
			if(twos>max_gc) max_gc=twos;
			if(twos<min_gc) min_gc=twos;
			count++;
		}
		ranges[0]=min_subgc;
		ranges[1]=max_subgc;
		ranges[2]=min_supergc;
		ranges[3]=max_supergc;
		ranges[4]=min_gc;
		ranges[5]=max_gc;
	
		return(count);
	}
}

static void go_to_index(const unsigned int index, const unsigned int dimensions, unsigned int *distribution, unsigned int *ranges){

	/* To verify success, confirm that distribution[0]!=NULL_INDEX */
	unsigned int count=0;
	
	unsigned int zeros=0;
	unsigned int ones=0;
	unsigned int twos=0;
	{
		int i;
		for(i=0; i<dimensions; i++){
			distribution[i]=NULL_INDEX;
		}
		while(increment_index(dimensions, distribution, ranges, &zeros, &ones, &twos)!=0){
			if(index==count) return;
			count++;
		}
	}		
}

void AIR_free_ica_sources_aic15(void ***constants_ptr){

	/* 
	 * This frees all the memory allocated by ica_init, ica_svd_init and ica_aic_logP_init
	 *
	 * If you don't want to free some variable, put a reference to it in some
	 * other variable and set the corresponding value in constants to NULL
	 * before calling this subroutine
	 */
	 
	void **constants=*constants_ptr;
	if(*constants_ptr==NULL) return;

	// constants[0] is pointer to ica_p
	// constants[1] is pointer to ica_n
	if(constants[2]) ica_pca_free_2((void *)constants[2]);	// ica_x
	if(constants[3]) ica_pca_free_2((void *)constants[3]);	// ica_w
	if(constants[4]) ica_pca_free_2((void *)constants[4]);	// ica_s
	if(constants[5]) ica_pca_free(constants[5]);					// ica_ipvt
	if(constants[6]) ica_pca_free_2((void *)constants[6]);	// ica_dw
	if(constants[7]) ica_pca_free(constants[7]);					// ica_work
	// constants[9] is pointer to subgc
	// constants[10] is pointer to supergc
	// ica_p-subgc-supergc is gc in subroutines that use constants
	// constants[11] is pointer to log of absolute value of determinant of W

	if(constants[12]) ica_pca_free_2((void *)constants[12]);	// ica_x_copy

	if(constants[14]) ica_pca_free_2((void *)constants[14]);	// qr_matrix
	if(constants[15]) ica_pca_free(constants[15]);					// qraux
	if(constants[16]) ica_pca_free(constants[16]);					// null_matrix_column
	if(constants[17]) ica_pca_free(constants[17]);					// ica_svd_work
	if(constants[18]) ica_pca_free_2((void *)constants[18]);	// s2
	if(constants[19]) ica_pca_free(constants[19]);					// variance
	
	// The following items are not used by ucminf calls
	if(constants[40]) ica_pca_free(constants[40]);					// internal ica_p and ica_n
	if(constants[41]) ica_pca_free(constants[41]);					// distribution
	if(constants[42]) ica_pca_free(constants[42]);					// probability
	if(constants[43]) ica_pca_free(constants[43]);					// distribution_aic_best
	if(constants[44]) ica_pca_free(constants[44]);					// kpvt
	if(constants[47]) ica_pca_free(constants[47]);					// ica_ica_w
	if(constants[48]) ica_pca_free(constants[48]);					// ica_ica_work
	if(constants[50]) ica_pca_free_2((void *)constants[50]);	// null_matrix
	
	// This item is needed for bias computation and is therefore just stored here
	if(constants[51]) ica_pca_free(constants[51]);					// ica_scale_x 
	
	// This item is needed for cross validation and is therefore just stored here
	if(constants[52]) ica_pca_free(constants[52]);
	
	// These are used for the initial svd
	if(constants[80]) ica_pca_free(constants[80]);					// internal n and points
	if(constants[81]) ica_pca_free(constants[81]);					// internal n and points for em_transpose
	if(constants[82]) ica_pca_free_2(constants[82]);			// em_transpose
	if(constants[83]) ica_pca_free(constants[83]);					// v
	if(constants[84]) ica_pca_free(constants[84]);					// work
	if(constants[85]) ica_pca_free(constants[85]);					// e
	
	
	// These are used for aic_logP
	if(constants[90]) ica_pca_free(constants[90]);					// internal enough_permutations
	if(constants[91]) ica_pca_free(constants[91]);					// aic_logP
	
	if(constants) ica_pca_free(constants);
	*constants_ptr=NULL;
}

static AIR_Error AIR_ica_svd_init(void ***constants, const unsigned int n, const unsigned int points){

	AIR_Boolean n_is_too_small=FALSE;
	AIR_Boolean points_is_too_small=FALSE;
	
	void **svd_constants;
	
	if(!*constants){
		// Allocate constants and initialize to all NULL
		*constants=ica_pca_malloc((size_t)CONSTANT_COUNT, sizeof(void **));
		if(!*constants){
			return AIR_MEMORY_ALLOCATION_ERROR;;
		}
		{
			unsigned int i;
			for(i=0; i<CONSTANT_COUNT; i++){
				(*constants)[i]=NULL;
			}
		}
	}
	svd_constants=*constants;
	
	{
		unsigned int *temp;
		
		if(!svd_constants[80]){
		
			svd_constants[80]=ica_pca_malloc((size_t)2,sizeof(unsigned int));
			if(!svd_constants[80]) return AIR_MEMORY_ALLOCATION_ERROR;
				
			temp=(unsigned int *)(svd_constants[80]);
			temp[0]=0;
			temp[1]=0;
		}
		else temp=(unsigned int *)(svd_constants[80]);
		
		n_is_too_small=(temp[0]<n);
		points_is_too_small=(temp[1]<points);
	}
	if(n>points){
		AIR_Boolean n2_is_too_small=FALSE;
		AIR_Boolean points2_is_too_small=FALSE;
		{
			unsigned int *temp;
			if(!svd_constants[81]){
			
				svd_constants[81]=ica_pca_malloc((size_t)2,sizeof(unsigned int));
				if(!svd_constants[81]) return AIR_MEMORY_ALLOCATION_ERROR;
				
				temp=(unsigned int *)(svd_constants[81]);
					
				temp[0]=0;
				temp[1]=0;
			}
			else temp=(unsigned int *)(svd_constants[81]);
			
			n2_is_too_small=(temp[0]<n);
			points2_is_too_small=(temp[1]<points);
		}
		if(n2_is_too_small || points2_is_too_small || svd_constants[82]==NULL){
			
			if(svd_constants[82]){
				ica_pca_free_2(svd_constants[82]);
				svd_constants[82]=NULL;
			}
			{
				double **em_transpose=ica_pca_matrix2(n, points);
				if(!em_transpose) return AIR_MEMORY_ALLOCATION_ERROR;
				
				svd_constants[82]=em_transpose;
			}
		}
		{
			unsigned int *temp=(unsigned int *)(svd_constants[81]);
			
			if(n2_is_too_small) temp[0]=n;
			if(points2_is_too_small) temp[1]=points;
		}
	}
	
	if(n_is_too_small || points_is_too_small || svd_constants[83]==NULL){
	
		if(svd_constants[83]){
			ica_pca_free(svd_constants[83]);
			svd_constants[83]=NULL;
		}
		{
			double *v;
			unsigned int temp=points+1;
			if(n+1<temp) temp=n+1;
			v=ica_pca_matrix1(temp);
			if(!v) return AIR_MEMORY_ALLOCATION_ERROR;
			
			svd_constants[83]=v;
		}
	}
	if(points_is_too_small || svd_constants[84]==NULL){
	
		if(svd_constants[84]){
			ica_pca_free(svd_constants[84]);
			svd_constants[84]=NULL;
		}
		{
			double *work=ica_pca_matrix1(points);
			if(!work) return AIR_MEMORY_ALLOCATION_ERROR;
			
			svd_constants[84]=work;
		}
	}
	if(n_is_too_small || svd_constants[85]==NULL){
	
		if(svd_constants[85]){
			ica_pca_free(svd_constants[85]);
			svd_constants[85]=NULL;
		}
		{
			double *e=ica_pca_matrix1(n);
			if(!e) return AIR_MEMORY_ALLOCATION_ERROR;
			
			svd_constants[85]=e;
		}
	}
	{
		unsigned int *temp=(unsigned int *)(svd_constants[80]);
		
		if(n_is_too_small) temp[0]=n;
		if(points_is_too_small) temp[1]=points;
	}
	return 0;
}

static AIR_Error AIR_ica_aic_logP_init(void ***constants, const unsigned int total_permutations){

	AIR_Boolean total_permutations_is_too_small=FALSE;
	
	void **options_constants;
	
	if(!*constants){
		// Allocate constants and initialize to all NULL
		*constants=ica_pca_malloc((size_t)CONSTANT_COUNT, sizeof(void **));
		if(!*constants){
			return AIR_MEMORY_ALLOCATION_ERROR;;
		}
		{
			unsigned int i;
			for(i=0; i<CONSTANT_COUNT; i++){
				(*constants)[i]=NULL;
			}
		}
	}

	options_constants=*constants;
	
	{
		unsigned int *temp;
	
		if(!options_constants[90]){
		
			options_constants[90]=ica_pca_malloc((size_t)1, sizeof(unsigned int));
			if(!options_constants[90]) return AIR_MEMORY_ALLOCATION_ERROR;
				
			temp=(unsigned int *)(options_constants[90]);
			temp[0]=0;
		}
		else temp=(unsigned int *)(options_constants[90]);
		
		total_permutations_is_too_small=(temp[0]<total_permutations);
	}
		
	if(total_permutations_is_too_small || options_constants[91]==NULL){
	
		if(options_constants[91]){
			ica_pca_free(options_constants[91]);
			options_constants[91]=NULL;
		}
		{
			struct hoods *aic_logP=(struct hoods *)ica_pca_malloc((size_t)total_permutations, sizeof(struct hoods));
			if(!aic_logP) return AIR_MEMORY_ALLOCATION_ERROR;
			
			options_constants[91]=aic_logP;
		}
	}
	{
		unsigned int *temp=(unsigned int *)(options_constants[90]);
		
		if(total_permutations_is_too_small) temp[0]=total_permutations;
	}
	return 0;
}

static AIR_Error AIR_ica_init(void ***constants, const unsigned int ica_p, 
	const unsigned int ica_n, const unsigned int hinted_sources, 
	AIR_Boolean xval_and_do_not_leave_uncentered, AIR_Boolean randomize_ica_x){

	/*
	 * This will allocate all the memory needed to perform ica_optimization
	 * using a call to ucminf
	 *
	 * If *constants is NULL, it will allocate it as well
	 *
	 * If it is possible that ica_p or ica_n have changed, this routine can
	 * be called again--it will free and reallocate any variables that don't 
	 * have sufficient storage
	 *
	 * *constants[40] is a pointer to an unsigned int array of length 2 that 
	 * stores ica_p and ica_n
	 *
	 * if *constants[40] is NULL, everything must be allocated, otherwise, it 
	 * is only necessary to reallocate when ica_p or ica_n have changed
	 *
	 * If an error is returned, it is safe to call ica_term or to re-call
	 * this routine
	 */
	 
	AIR_Boolean ica_p_is_too_small=FALSE;
	AIR_Boolean ica_n_is_too_small=FALSE;
	
	void **ica_constants;
	 
	 
	if(!*constants){
		// Allocate constants and initialize to all NULL
		*constants=ica_pca_malloc(CONSTANT_COUNT, sizeof(void **));
		if(!*constants){
			return AIR_MEMORY_ALLOCATION_ERROR;;
		}
		{
			unsigned int i;
			for(i=0; i<CONSTANT_COUNT; i++){
				(*constants)[i]=NULL;
			}
		}
	}

	ica_constants=*constants;

	{
		unsigned int *temp;
	
		if(!ica_constants[40]){
		
			ica_constants[40]=ica_pca_malloc((size_t)2, sizeof(unsigned int));
			if(!ica_constants[40]) return AIR_MEMORY_ALLOCATION_ERROR;
			
			temp=(unsigned int *)(ica_constants[40]);
			temp[0]=0;
			temp[1]=0;
		}
		else temp=(unsigned int *)(ica_constants[40]);
		
		ica_p_is_too_small=(temp[0]<ica_p);
		ica_n_is_too_small=(temp[1]<ica_n);
	}
	
	// ica_constants[0] is address of ica_p;
	// ica_constants[1] is address of ica_n;
	
	if(ica_p_is_too_small || ica_n_is_too_small || ica_constants[2]==NULL){
	
		if(ica_constants[2]){
			ica_pca_free_2(ica_constants[2]);
			ica_constants[2]=NULL;
		}
		{
			double **ica_x=ica_pca_matrix2(ica_p, ica_n);
			if(!ica_x) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[2]=ica_x;
		}
	}
	
	if(ica_p_is_too_small || ica_constants[3]==NULL){
	
		if(ica_constants[3]){
			ica_pca_free_2(ica_constants[3]);
			ica_constants[3]=NULL;
		}
		{
			double **ica_w=ica_pca_matrix2(ica_p, ica_p);
			if(!ica_w) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[3]=ica_w;
		}
	}
	
	if(ica_p_is_too_small || ica_n_is_too_small || ica_constants[4]==NULL){

		if(ica_constants[4]){
			ica_pca_free_2(ica_constants[4]);
			ica_constants[4]=NULL;
		}
		{
			double **ica_s=ica_pca_matrix2(ica_p, ica_n);
			if(!ica_s) return AIR_MEMORY_ALLOCATION_ERROR;
	
			ica_constants[4]=ica_s;
		}
	}
	
	if(ica_p_is_too_small || ica_constants[5]==NULL){
	
		if(ica_constants[5]){
			ica_pca_free(ica_constants[5]);
			ica_constants[5]=NULL;
		}
		{
			unsigned int *ica_ipvt=(unsigned int *)ica_pca_malloc((size_t)ica_p, sizeof(unsigned int));
			if(!ica_ipvt) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[5]=ica_ipvt;
		}
	}
	
	if(ica_p_is_too_small || ica_constants[6]==NULL){

		if(ica_constants[6]){
			ica_pca_free_2(ica_constants[6]);
			ica_constants[6]=NULL;
		}
		{
			double **ica_dw=ica_pca_matrix2(ica_p, ica_p);
			if(!ica_dw) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[6]=ica_dw;
		}
	}
	
	if(ica_p_is_too_small || ica_constants[7]==NULL){

		if(ica_constants[7]){
			ica_pca_free(ica_constants[7]);
			ica_constants[7]=NULL;
		}
		{
			double *ica_work=(double *)ica_pca_malloc((size_t)ica_p, sizeof(double));
			if(!ica_work) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[7]=ica_work;
		}
	}
	
	// ica_constants[9]=subgc;
	// ica_constants[10]=supergc;
	// gc=m-subgc-supergc
	
	if(xval_and_do_not_leave_uncentered || randomize_ica_x){
		if(ica_p_is_too_small || ica_n_is_too_small || ica_constants[12]==NULL){
		
			if(ica_constants[12]){
				ica_pca_free_2(ica_constants[12]);
				ica_constants[12]=NULL;
			}
			{
				double **ica_x_copy=ica_pca_matrix2(ica_p, ica_n);
				if(!ica_x_copy) return AIR_MEMORY_ALLOCATION_ERROR;
				
				ica_constants[12]=ica_x_copy;
			}
		}
	}

	
	if(ica_p_is_too_small || ica_constants[14]==NULL){
	
		if(ica_constants[14]){
			ica_pca_free_2(ica_constants[14]);
			ica_constants[14]=NULL;
		}
		{
			double **qr_matrix=ica_pca_matrix2(ica_p, ica_p);
			if(!qr_matrix) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[14]=qr_matrix;
		}
	}
		
	
	if(ica_p_is_too_small || ica_constants[15]==NULL){
			
		if(ica_constants[15]){
			ica_pca_free(ica_constants[15]);
			ica_constants[15]=NULL;
		}
		{
			double *qraux=(double *)ica_pca_malloc((size_t)ica_p, sizeof(double));
			if(!qraux) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[15]=qraux;
		}
	}
	
	if(ica_n_is_too_small || ica_constants[16]==NULL){
	
		if(ica_constants[16]){
			ica_pca_free(ica_constants[16]);
			ica_constants[16]=NULL;
		}
		{
			double *null_matrix_column=(double *)ica_pca_malloc((size_t)(ica_n+1), sizeof(double));
			if(!null_matrix_column) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[16]=null_matrix_column;
		}
	}
	if(ica_n_is_too_small || ica_constants[17]==NULL){
		
		if(ica_constants[17]){
			ica_pca_free(ica_constants[17]);
			ica_constants[17]=NULL;
		}
		{
			double *ica_svd_work=(double *)ica_pca_malloc((size_t)ica_n, sizeof(double));
			if(!ica_svd_work) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[17]=ica_svd_work;
		}
	}
	
	if(ica_p_is_too_small || ica_constants[18]==NULL){
		
		if(ica_constants[18]){
			ica_pca_free_2(ica_constants[18]);
			ica_constants[18]=NULL;
		}
		{
			double **ica_s2=ica_pca_matrix2(ica_n, ica_p);
			if(!ica_s2) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[18]=ica_s2;
		}
	}

	if(ica_p_is_too_small || ica_constants[19]==NULL){
		
		if(ica_constants[19]){
			ica_pca_free_2((void *)ica_constants[19]);
			ica_constants[19]=NULL;
		}
		{
			double *variance=(double *)ica_pca_malloc((size_t)ica_p, sizeof(double));
			if(!variance) return AIR_MEMORY_ALLOCATION_ERROR;

			ica_constants[19]=variance;
		}
	}
	
	if(ica_p_is_too_small || ica_constants[41]==NULL){
	
		if(ica_constants[41]){
			ica_pca_free(ica_constants[41]);
			ica_constants[41]=NULL;
		}
		{
			unsigned int *distribution=(unsigned int *)ica_pca_malloc((size_t)ica_p, sizeof(unsigned int));
			if(!distribution) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[41]=distribution;
		}
	}
	if(ica_p_is_too_small || ica_constants[42]==NULL){
		
		if(ica_constants[42]){
			ica_pca_free(ica_constants[42]);
			ica_constants[42]=NULL;
		}
		{
			double *probability=(double *)ica_pca_malloc((size_t)(ica_p+1), sizeof(double));
			if(!probability) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[42]=probability;
		}
	}
	if(ica_p_is_too_small || ica_constants[43]==NULL){
	
		if(ica_constants[43]){
			ica_pca_free(ica_constants[43]);
			ica_constants[43]=NULL;
		}
		{
			unsigned int *distribution_aic_best=(unsigned int *)ica_pca_malloc((size_t)ica_p, sizeof(unsigned int));
			if(!distribution_aic_best) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[43]=distribution_aic_best;
		}
	}
	if(ica_p_is_too_small || ica_constants[44]==NULL){
	
		if(ica_constants[44]){
			ica_pca_free(ica_constants[44]);
			ica_constants[44]=NULL;
		}
		{
			unsigned int *kpvt=(unsigned int *)ica_pca_malloc((size_t)ica_p, sizeof(unsigned int));
			if(!kpvt) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[44]=kpvt;
		}
	}

	if(ica_p_is_too_small || ica_constants[47]==NULL){
	
		if(ica_constants[47]){
			ica_pca_free(ica_constants[47]);
			ica_constants[47]=NULL;
		}
		{
			double *ica_ica_w=(double *)ica_pca_malloc((size_t)ica_p*ica_p, sizeof(double));
			if(!ica_ica_w) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[47]=ica_ica_w;
		}
	}
	if(ica_p_is_too_small || ica_constants[48]==NULL){
	
		if(ica_constants[48]){
			ica_pca_free(ica_constants[48]);
			ica_constants[48]=NULL;
		}
		{
			double *ica_ica_work=(double *)ica_pca_malloc(((size_t)ica_p*ica_p*(ica_p*ica_p+13)/2), sizeof(double));
			if(!ica_ica_work) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[48]=ica_ica_work;
		}
	}
	if(hinted_sources){
		if(ica_p_is_too_small || ica_constants[50]==NULL){
		
			if(ica_constants[50]){
				ica_pca_free_2((void *)ica_constants[50]);
				ica_constants[50]=NULL;
			}
			{
				double **null_matrix=ica_pca_matrix2(ica_p, ica_p);
				if(!null_matrix) return AIR_MEMORY_ALLOCATION_ERROR;
				
				ica_constants[50]=null_matrix;
			}
		}
	}
	if(ica_constants[51]==NULL){
		
		double *ica_scale_x=ica_pca_matrix1(sizeof(double));
		if(!ica_scale_x) return AIR_MEMORY_ALLOCATION_ERROR;
		
		ica_constants[51]=ica_scale_x;
	}
	if(ica_p_is_too_small || ica_constants[47]==NULL){
	
		if(ica_constants[52]){
			ica_pca_free(ica_constants[52]);
			ica_constants[52]=NULL;
		}
		{
			double *ica_ica_w2=(double *)ica_pca_malloc((size_t)ica_p*ica_p, sizeof(double));
			if(!ica_ica_w2) return AIR_MEMORY_ALLOCATION_ERROR;
			
			ica_constants[52]=ica_ica_w2;
		}
	}

	{
		unsigned int *temp=(unsigned int *)(ica_constants[40]);
		
		if(ica_p_is_too_small) temp[0]=ica_p;
		if(ica_n_is_too_small) temp[1]=ica_n;
	}
	
	return 0;
}


AIR_Error AIR_ica_sources_aic15(const unsigned int seed, const unsigned int offset_random, unsigned int *sample_order, const unsigned int sample_offset, const unsigned int samples, void ***constants_ptr, const unsigned int inf_crit, const unsigned int fold, const double xval_epsilon, unsigned int *ranges, double **em, const unsigned int points, const unsigned int n, const unsigned int hinted_subgaussian_sources, const unsigned int hinted_supergaussian_sources, const unsigned int hinted_unspecified_sources, double **source_matrix, unsigned int *components, double ***ica_s_ptr, unsigned int desired_initialization, unsigned int **distribution_ptr, double **variance_ptr, double **probability_ptr, double *loglikelihood, const char *outprefix, AIR_Error (*data_saver)(const AIR_Boolean, const char *, const unsigned int, const unsigned int, const unsigned int, double **, double, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const unsigned int, double *, void **, const AIR_Boolean), void **plot_constants, int *model_of_interest, double *moi_probability, const AIR_Boolean leave_rows_uncentered, const AIR_Boolean verbose, const AIR_Boolean ow){
 	
 	void **constants;
 	 	
	*ica_s_ptr=NULL;
	*distribution_ptr=NULL;
	*variance_ptr=NULL;
	*probability_ptr=NULL;
	
	unsigned int stride=1;
	
	const unsigned int hinted_sources=hinted_subgaussian_sources+hinted_supergaussian_sources+hinted_unspecified_sources;
	const unsigned int nonpermuted_sources=hinted_subgaussian_sources+hinted_supergaussian_sources;
	
#ifndef USING_R
	if(seed!=0){
		srand(seed);
	}
#endif
		
	if(points==0 || n==0 || *components==0) return 0;	// Nothing to do
	
	if(inf_crit==3 && fold){
		if(fold==1){
#ifndef USING_R
			printf("%s: %d: ",__FILE__,__LINE__);
			printf("cross-validation fold of 1 is not valid\n");
#else
			REprintf("%s: %d: ",__FILE__,__LINE__);
			REprintf("cross_validation fold of 1 is not valid\n");
#endif
			return AIR_USER_INTERFACE_ERROR;
		}
		if(samples!=0 && samples!=n){
#ifndef USING_R
			printf("%s: %d: ",__FILE__,__LINE__);
			printf("all %u samples must be used when %u fold cross-validation is specified\n", n, fold);
#else
			REprintf("%s: %d: ",__FILE__,__LINE__);
			REprintf("all %u samples must be used when %u fold cross-validation is specified\n", n, fold);
#endif
			return AIR_USER_INTERFACE_ERROR;
		}
		stride=n/fold;
		if(stride*fold!=n){
#ifndef USING_R
			printf("%s: %d: ",__FILE__,__LINE__);
			printf("cross-validation fold (%u) is not a factor of n (%u)\n", fold, n);
#else
			REprintf("%s: %d: ",__FILE__,__LINE__);
			REprintf("cross_validation fold (%u) is not a factor of n (%u)\n", fold, n);
#endif
			return AIR_USER_INTERFACE_ERROR;
		}
	}
	
	// Look for fatal errors that should have been caught by user interface
	if( (hinted_sources>0) && source_matrix==NULL){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("Programming error: source_matrix cannot be NULL when hinted_sources are non-zero\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("Programming error: source_matrix cannot be NULL when hinted_sources are non-zero\n");
#endif
		return AIR_USER_INTERFACE_ERROR;
	}
	if(inf_crit>3 && leave_rows_uncentered){
#ifndef USING_R
		printf("%s: %d: ",__FILE__,__LINE__);
		printf("Cannot leave rows uncentered when using second order corrected AIC for gaussian sources\n");
#else
		REprintf("%s: %d: ",__FILE__,__LINE__);
		REprintf("Cannot leave rows uncentered when using second order corrected AIC for gaussian sources\n");
#endif
		return AIR_USER_INTERFACE_ERROR;
	}

	if(inf_crit==3 && samples){
		if(sample_offset+samples>n){
#ifndef USING_R
			printf("%s: %d: ",__FILE__,__LINE__);
			printf("Sum of sample offsets and number of samples cannot exceed number of observations\n");
#else
			REprintf("%s: %d: ",__FILE__,__LINE__);
			REprintf("Sum of sample offsets and number of samples cannot exceed number of observations\n");
#endif
			return AIR_USER_INTERFACE_ERROR;
		}
	}
	if(inf_crit==3 && sample_order){
		{
			unsigned int j;
			for(j=n; j>0; j--){
				int found=0;
				unsigned int i;
				for(i=0; i<n; i++){
					if(sample_order[i]==j){
						found=1;
						break;
					}
				}
				if(!found){
#ifndef USING_R
					printf("%s: %d: ",__FILE__,__LINE__);
					printf("Cross-validation sample order should contain the integers 1 through %u, the value %u is missing\n",n, j);
#else
					REprintf("%s: %d: ",__FILE__,__LINE__);
					REprintf("Cross-validation sample order should contain the integers 1 through %u, the value %u is missing\n",n, j);
#endif
					return AIR_USER_INTERFACE_ERROR;
				}
			}
		}
	}
	{
		AIR_Boolean use_subgc=(AIR_Boolean)TRUE;
		AIR_Boolean use_supergc=(AIR_Boolean)TRUE;
		AIR_Boolean use_gc=(AIR_Boolean)TRUE;
		
		if(ranges[1]<ranges[0] || ranges[1]==0) use_subgc=(AIR_Boolean)FALSE;
		if(ranges[3]<ranges[2] || ranges[3]==0) use_supergc=(AIR_Boolean)FALSE;
		if(ranges[5]<ranges[4] || ranges[5]==0) use_gc=(AIR_Boolean)FALSE;
		
		if(!use_subgc && !use_supergc && !use_gc){
#ifndef USING_R
			printf("%s: %d: ",__FILE__,__LINE__);
			printf("No models fit specified range restrictions\n");
#else
			REprintf("%s: %d: ",__FILE__,__LINE__);
			REprintf("No models fit specified range restrictions\n");
#endif
			return 0; // Nothing to do
		}
				
		if(outprefix!=NULL && data_saver!=NULL){
			// Run pre-analysis check of file saving
			AIR_Error errcode=data_saver((AIR_Boolean)TRUE, outprefix, 0, n, *components, NULL, 0, 0, 0, 0, 0, 0, 0, inf_crit, NULL, plot_constants, ow);
			if(errcode!=0) return errcode;
		}
		{
			if(!leave_rows_uncentered){
				/* Center each row by subtracting the row mean 
				 *
				 * The rows produced by SVD will also end up centered
				 *  but doing the centering before SVD assures
				 *  that the number components ultimately
				 *  computed will take the redundancy among row.
				 *  in account. For example, a square matrix
				 *  of size N will only allow N-1 sources
				 *  after centering.
				 */
				unsigned int i;
				
				for(i=0; i<points; i++){
				
					double row_mean=0.0;
					unsigned int j;
					
					for(j=0; j<n; j++){
				
						row_mean+=em[j][i];
					}
					row_mean/=n;
					for(j=0;j<n; j++){
					
						em[j][i]-=row_mean;
					}
					/* Re-average to make values sum as closely to zero as possible */
					row_mean=0;
					for(j=0; j<n; j++){
				
						row_mean+=em[j][i];
					}
					row_mean/=n;
					for(j=0;j<n; j++){
					
						em[j][i]-=row_mean;
					}
				}
			}
			{
				AIR_Error errcode=AIR_ica_svd_init(constants_ptr, n, points);
				if(errcode){
					return errcode;
				}
				constants=*constants_ptr;
			}
			{
				double *v=NULL;
				double **d=em;
				AIR_Boolean transpose=(AIR_Boolean)FALSE;
				double **em_transpose=NULL;
				
				if(n>points) transpose=(AIR_Boolean)TRUE;
				
				if(transpose){
					em_transpose=constants[82];
					{
						unsigned int j;
						
						for(j=0;j<points;j++){
						
							unsigned int i;
							
							for(i=0;i<n;i++){
							
								em_transpose[j][i]=em[i][j];
							}
						}
					}
				}
				/* Perform SVD */
				{
					v=constants[83];
									
					{
						double *work=constants[84];
						double *e=constants[85];
						{
							AIR_Error errcode;
							
							if(transpose){
								errcode=AIR_dsvdc(em_transpose,n,points,v,work,em_transpose,NULL,e,(unsigned int)2, (AIR_Boolean)FALSE,NULL);
							}
							else{
								errcode=AIR_dsvdc(em,points,n,v,e,NULL,d,work,(unsigned int)0,(AIR_Boolean)1,NULL);
							}
							if(errcode!=0){
#ifndef USING_R
								printf("%s: %d: ",__FILE__,__LINE__);
								printf("Error occurred in computing singular value decomposition\n");
#else
								REprintf("%s: %d: ",__FILE__,__LINE__);
								REprintf("Error occurred in computing singular value decomposition\n");
#endif
								AIR_free_ica_sources_aic15(constants_ptr);
								return(errcode);
							}
						}
						if(transpose){
							{
								unsigned int j;
								
								for(j=0;j<points;j++){
									
									unsigned int i;
									
									for(i=0;i<n;i++){
										// Note that d==em
										d[i][j]=em_transpose[j][i];
									}
								}
							}
						}
						/* Perform ICA */
						{
							unsigned int ica_n=n;	// Number of columns
							unsigned int ica_p=n;	// min(rows,columns-1,nontrivial components,requested components)
							double total_variance=0.0;
							
							if(points<ica_p) ica_p=points;
							
							/* Sum total variance */
							{
								int i;
								
								for(i=0; i<points && i<n; i++){
									
									double added_variance=v[i]*v[i];
									total_variance+=added_variance;
								}
							}
							/* Check for trivial components, decrementing ica_p to the number of non-trivial ones */
							{
								unsigned int first_trivial_index;
								double cumulated_variance=v[0]*v[0];
								
								for(first_trivial_index=1;first_trivial_index<ica_p;first_trivial_index++){
									{
										double added_variance=v[first_trivial_index]*v[first_trivial_index];
										if(added_variance<2.0*cumulated_variance*DBL_EPSILON) break;
										cumulated_variance+=added_variance;
									}
								}
								ica_p=first_trivial_index;
							}
							if(*components<ica_p) ica_p=*components;
							/* Cross validation runs decrement number of columns, so decrement ica_p so rows<=columns-1 */
							if(inf_crit==3){
								if(ica_p>ica_n-stride) ica_p=ica_n-stride;
							}
							
							*components=ica_p;
														
							{
								// Note that summing the components alone is a problem
								// if the upper range is UINT_MAX, could overflow
								// back to small numbers, hence check the individual
								// elements and the sum
								AIR_Boolean adequate_components=(AIR_Boolean)FALSE;
								unsigned int available=0;
								if(use_subgc){
									if(ranges[1]>=ica_p) adequate_components=(AIR_Boolean)TRUE;
									available+=ranges[1];
								}
								if(use_supergc){
									if(ranges[3]>=ica_p) adequate_components=(AIR_Boolean)TRUE;
									available+=ranges[3];
								}
								if(use_gc){
									if(ranges[5]>=ica_p) adequate_components=(AIR_Boolean)TRUE;
									available+=ranges[5];
								}
								if(available>=ica_p) adequate_components=(AIR_Boolean)TRUE;
								if(!adequate_components){
#ifndef USING_R
									printf("%s: %d: ",__FILE__,__LINE__);
									printf("specified ranges are too low for the number of components\n");
#else
									REprintf("%s: %d: ",__FILE__,__LINE__);
									REprintf("specified ranges are too low for the number of components\n");
#endif
									AIR_free_ica_sources_aic15(constants_ptr);
									return AIR_USER_INTERFACE_ERROR;
								}
							}
							{
								unsigned int mandatory=0;
								if(use_subgc){
									mandatory+=ranges[0];
								}
								if(use_supergc){
									mandatory+=ranges[2];
								}
								if(use_gc){
									mandatory+=ranges[4];
								}
								if(mandatory>ica_p){
#ifndef USING_R
									printf("%s: %d: ",__FILE__,__LINE__);
									printf("specified ranges are too high for the number of components\n");
#else
									REprintf("%s: %d: ",__FILE__,__LINE__);
									REprintf("specified ranges are too high for the number of components\n");
#endif
									AIR_free_ica_sources_aic15(constants_ptr);
									return AIR_USER_INTERFACE_ERROR;
								}
							}
							if(hinted_sources>ica_p){
#ifndef USING_R
								printf("%s: %d: ",__FILE__,__LINE__);
								printf("Number of initialized sources (%u) cannot exceed the number of retained components (%u)\n", hinted_sources, ica_p);
#else
								REprintf("%s: %d: ",__FILE__,__LINE__);
								REprintf("Number of initialized sources (%u) cannot exceed the number of retained components (%u)\n", hinted_sources, ica_p);
#endif
								AIR_free_ica_sources_aic15(constants_ptr);
								return AIR_ICA_EXCESS_SOURCE_ERROR;
							}
							if(hinted_subgaussian_sources>ranges[1]){
#ifndef USING_R
								printf("%s: %d: ",__FILE__,__LINE__);
								printf("Number of hinted subgaussian sources (%u) exceeds maximum number of subgaussian sources (%u)\n", hinted_subgaussian_sources, ranges[1]);
#else
								REprintf("%s: %d: ",__FILE__,__LINE__);
								REprintf("Number of hinted subgaussian sources (%u) exceeds maximum number of subgaussian sources (%u)\n", hinted_subgaussian_sources, ranges[1]);
#endif
								AIR_free_ica_sources_aic15(constants_ptr);
								return AIR_ICA_EXCESS_SOURCE_ERROR;
							}
							if(hinted_supergaussian_sources>ranges[3]){
#ifndef USING_R
								printf("%s: %d: ",__FILE__,__LINE__);
								printf("Number of hinted supergaussian sources (%u) exceeds maximum number of supergaussian sources (%u)\n", hinted_supergaussian_sources, ranges[3]);
#else
								REprintf("%s: %d: ",__FILE__,__LINE__);
								REprintf("Number of hinted supergaussian sources (%u) exceeds maximum number of supergaussian sources (%u)\n", hinted_supergaussian_sources, ranges[3]);
#endif
								AIR_free_ica_sources_aic15(constants_ptr);
								return AIR_ICA_EXCESS_SOURCE_ERROR;
							}
							{
								AIR_Error errcode=AIR_ica_init(constants_ptr, ica_p, ica_n, hinted_sources, (sample_order || !leave_rows_uncentered || fold) && (inf_crit==3), (seed!=0 && hinted_sources==0));
								if(errcode){
									AIR_free_ica_sources_aic15(constants_ptr);
									return errcode;
								}
								constants=*constants_ptr;
							}
							{
								double **ica_x=constants[2];
								//double **ica_w=constants[3];
								double **ica_s=constants[4];

								double **qr_matrix=constants[14];
								
								double *variance=constants[19];
								
								unsigned int *distribution=constants[41];
								double *probability=constants[42];
								unsigned int *distribution_aic_best=constants[43];
								unsigned int *kpvt=constants[44];
								
								double *ica_ica_w=constants[47];
								double *ica_ica_work=constants[48];
								
								double **null_matrix=constants[50];
								
								const unsigned int all_sources=hinted_sources;
								
								unsigned int subgc=0;
								unsigned int supergc=0;
								unsigned int gc=0;
								
								double detw=0.0;
														
								{
									double ica_scale_x=0.0;
									unsigned int permuted_ranges[6];
									{
										/* Adjust ranges to account for hinted sources of specified distribution */
										{
											unsigned int i=hinted_subgaussian_sources;
											permuted_ranges[0]=ranges[0];
											while(i>0 && permuted_ranges[0]>0){
												i--;
												permuted_ranges[0]--;
											}
										}
										{
											unsigned int i=hinted_subgaussian_sources;
											permuted_ranges[1]=ranges[1];
											while(i>0 && permuted_ranges[1]>0){
												i--;
												permuted_ranges[1]--;
											}
										}
										{
											unsigned int i=hinted_supergaussian_sources;
											permuted_ranges[2]=ranges[2];
											while(i>0 && permuted_ranges[2]>0){
												i--;
												permuted_ranges[2]--;
											}
										}
										{
											unsigned int i=hinted_supergaussian_sources;
											permuted_ranges[3]=ranges[3];
											while(i>0 && permuted_ranges[3]>0){
												i--;
												permuted_ranges[3]--;
											}
										}
										permuted_ranges[4]=ranges[4];
										permuted_ranges[5]=ranges[5];
									}									
									constants[0]=&ica_p;
									constants[1]=&ica_n;
									constants[9]=&subgc;
									constants[10]=&supergc;
									/* constants[11]=&gc; */
									/* gc=ica_p-subgc-supergc */
									constants[11]=&detw;
																		
									/* Transpose and copy non-trivial scores into ica_x */
									if(transpose){
										unsigned int j;
										for(j=0;j<ica_n;j++){
											unsigned int i;
											for(i=0;i<ica_p;i++){
												ica_x[j][i]=v[i]*d[j][i];
												if(fabs(ica_x[j][i])>ica_scale_x) ica_scale_x=fabs(ica_x[j][i]);
											}
										}
									}
									else{
										unsigned int j;
										for(j=0;j<ica_n;j++){
											unsigned int i;
											for(i=0;i<ica_p;i++){
												ica_x[j][i]=v[i]*d[i][j];
												if(fabs(ica_x[j][i])>ica_scale_x) ica_scale_x=fabs(ica_x[j][i]);
											}
										}
									}
									if(all_sources!=0){
										{
											unsigned int j;
											for(j=0; j<all_sources; j++){
												
												unsigned int i;
												for(i=0; i<ica_n; i++){
													source_matrix[j][i]*=ica_scale_x;
												}
											}
										}

										/* Convert the sources into W' in the qr_matrix */
										{
											unsigned int j;
											for(j=0;j<all_sources;j++){
											
												unsigned int i;
												for(i=0;i<ica_p;i++){
												
													double sum=0.0;
													
													if(transpose){
													
														unsigned int k;
														
														for(k=0;k<ica_n;k++){
														
															sum+=source_matrix[j][k]*d[k][i];
														}			
													}
													else{
													
														unsigned int k;
														
														for(k=0;k<ica_n;k++){
														
															sum+=source_matrix[j][k]*d[i][k];
														}
													}
													qr_matrix[j][i]=sum/v[i];
												}
											}
										}
										
										// Find the null space
										{
											AIR_Error errcode=AIR_nullbasis(qr_matrix, ica_p, all_sources, null_matrix, NULL, NULL, NULL);
											if(errcode!=0){
												AIR_free_ica_sources_aic15(constants_ptr);
												return errcode;
											}
										}
										
										// Compute the sources into null_matrix
										{
											unsigned int j;
											for(j=0;j<all_sources;j++){
											
												unsigned int i;
												for(i=0;i<ica_p;i++){
												
													double sum=0.0;
													
													if(transpose){
													
														unsigned int k;
														
														for(k=0;k<ica_n;k++){
														
															sum+=source_matrix[j][k]*d[k][i];
														}	
													}
													else{
													
														unsigned int k;
														
														for(k=0;k<ica_n;k++){
														
															sum+=source_matrix[j][k]*d[i][k];
														}		
													}
													null_matrix[j][i]=sum/v[i];
												}
											}
										}				
									}
									
									/* Rescale ica_x */
									{
										unsigned int j;
										
										for(j=0;j<ica_n;j++){
										
											unsigned int i;
											
											for(i=0;i<ica_p;i++){
												ica_x[j][i]/=ica_scale_x;
											}
										}
									}
									if(seed!=0 && all_sources==0){
										/* Randomly rotate ica_x */
										double **r=AIR_matrix2(ica_p, ica_p);
										double ***wrk=AIR_matrix3(ica_p, ica_p, 4);
										double **ica_x_copy=constants[12];

										
										AIR_random_rotation(offset_random, r, ica_p, ica_n, ica_x, ica_x_copy, wrk);
										{
											unsigned int j;
											
											for(j=0; j<ica_n; j++){
											
												unsigned int i;
												
												for(i=0; i<ica_p; i++){
												
													ica_x[j][i]=ica_x_copy[j][i];
												}
											}
										}
										AIR_free_2(r);
										AIR_free_3(wrk);
									}
									else if( (sample_order || !leave_rows_uncentered || fold) && (inf_crit==3)){
									
										/* Make a backup copy of ica_x */
										double **ica_x_copy=constants[12];
										unsigned int j;
										
										for(j=0; j<ica_n; j++){
										
											unsigned int i;
											
											for(i=0; i<ica_p; i++){
											
												ica_x_copy[j][i]=ica_x[j][i];
											}
										}
									}
									if(seed!=0 && all_sources!=0){
										/* Randomly rotate null_matrix */
										double **r=AIR_matrix2(ica_p-all_sources, ica_p-all_sources);
										double ***wrk=AIR_matrix3(ica_p-all_sources, ica_p-all_sources, 4);
										double **tmp=AIR_matrix2(ica_p, ica_p-all_sources);
										
										AIR_random_rotation2(offset_random, r, ica_p, ica_p-all_sources, null_matrix+all_sources, tmp, wrk);
										AIR_free_2(r);
										AIR_free_3(wrk);
										AIR_free_2(tmp);
									}
									/* Store ica_scale_x in case it's needed for bias computations */
									{
										double *temp=(double *)constants[51];
										*temp=ica_scale_x;
									}
									/* Perform ICA optimization */
									{
										const unsigned int total_permutations=count_valid_indices(ica_p-nonpermuted_sources, distribution+nonpermuted_sources, permuted_ranges);
										unsigned int enough_permutations=total_permutations;
										unsigned int working_permutations=total_permutations;
										unsigned int min_possible_gaussians=0;
										unsigned int max_possible_gaussians=0;
										
										if(desired_initialization!=0 && desired_initialization<enough_permutations){
											enough_permutations=desired_initialization;
											working_permutations=desired_initialization;
										}
										
										if(inf_crit>3 && ranges[5]+3>=ica_n){
#ifndef USING_R
											printf("%s: %d: ",__FILE__,__LINE__);
											printf("When using second order corrected AIC for gaussian sources, the number of gaussian components plus three must be less than the number of columns\n");
#else
											REprintf("%s: %d: ",__FILE__,__LINE__);
											REprintf("When using second order corrected AIC for gaussian sources, the number of gaussian components plus three must be less than the number of columns\n");
#endif
											AIR_free_ica_sources_aic15(constants_ptr);
											return AIR_USER_INTERFACE_ERROR;
										}
										if(inf_crit>3 && (ranges[1]>0 || ranges[3]>0)){
#ifndef USING_R
											printf("%s: %d: ",__FILE__,__LINE__);
											printf("Second order corrected AIC can only be used when all sources are Gaussian\n");
#else
											REprintf("%s: %d: ",__FILE__,__LINE__);
											REprintf("Second order corrected AIC can only be used when all sources are Gaussian\n");
#endif
											AIR_free_ica_sources_aic15(constants_ptr);
											return AIR_USER_INTERFACE_ERROR;
										}
										{
											if(use_gc){
												max_possible_gaussians=ica_p;
											}
											if(!use_subgc && !use_supergc){
												min_possible_gaussians=ica_p;
											}
										}
										{
											AIR_Error errcode=AIR_ica_aic_logP_init(constants_ptr, enough_permutations);
											if(errcode){
												AIR_free_ica_sources_aic15(constants_ptr);
												return errcode;
											}
											constants=*constants_ptr;
										}
										{
											struct hoods *aic_logP=(struct hoods *)constants[91];
											if(desired_initialization==0){
												unsigned int index;
												for(index=0; index<enough_permutations; index++){
													aic_logP[index].rank=index;
												}
											}
											else{

												if(desired_initialization>total_permutations){
#ifndef USING_R
													printf("%s: %d: ",__FILE__,__LINE__);
													printf("desired initialization must be less than total permutations (%u)\n",total_permutations);
#else
													REprintf("%s: %d: ",__FILE__,__LINE__);
													REprintf("desired initialization must be less than total permutations (%u)\n",total_permutations);
#endif
													AIR_free_ica_sources_aic15(constants_ptr);
													return(AIR_USER_INTERFACE_ERROR);
												}
												if(verbose){
#ifndef USING_R
													printf("Restricting analysis to initialization %u/%u\n", desired_initialization, total_permutations);
#else
													Rprintf("Restricting analysis to initialization %u/%u\n", desired_initialization, total_permutations);
#endif
												}
												aic_logP[0].rank=desired_initialization-1;
												{
													unsigned int index;
													for(index=1; index<enough_permutations; index++){
														aic_logP[index].rank=enough_permutations;
														aic_logP[index].value=0.0;
													}
												}
											}
											{
												unsigned int pass;
												unsigned int passes=1;
												unsigned int working_inf_crit=inf_crit;
												
												if(inf_crit==3 && xval_epsilon!=0.0){
													passes=2;
													working_inf_crit=0;
												}
												
												for(pass=0; pass<passes; pass++){
												
													struct hoods *current_aic_logP=aic_logP;
													double aic_lowest=DBL_MAX;
													
													if(pass>0) working_inf_crit=3;
													{
														unsigned int index=0;
														unsigned int zeros=0;
														unsigned int ones=0;
														unsigned int twos=0;
														{
															unsigned int i=0;
															unsigned int i2;

															/* Mark hinted subgaussian sources */
															for(i2=0; i2<hinted_subgaussian_sources; i++, i2++){
																distribution[i]=0;
															}
															/* Mark hinted supergaussian sources */
															for(i2=0; i2<hinted_supergaussian_sources; i++, i2++){
																distribution[i]=1;
															}
															/* Initialize unknown sources */
															for(; i<ica_p; i++){
																distribution[i]=NULL_INDEX;
															}
														}
														if(verbose && !desired_initialization){
#ifndef USING_R
															if(pass>0){
																printf("testing %u distinct initializations of %u total initializations\n", working_permutations, total_permutations);
															}
															else if(passes>1){
																printf("searching %u initializations for distinct initializations\n", working_permutations);
															}
															else{
																printf("testing %u initializations of ICA algorithm\n", working_permutations);
															}
#else
															if(pass>0){
																Rprintf("testing %u distinct initializations of %u total initializations\n", working_permutations, total_permutations);
															}
															else if(passes>1){
																Rprintf("searching %u initializations for distinct initializations\n", working_permutations);
															}
															else{
																Rprintf("testing %u initializations of ICA algorithm\n", working_permutations);
															}
#endif
														}
														for(index=0; index<enough_permutations; index++){
																												
															if(current_aic_logP->rank==enough_permutations) break;
														
															/* Get the next distribution */
															increment_index(ica_p-nonpermuted_sources, distribution+nonpermuted_sources, permuted_ranges, &zeros, &ones, &twos);
															
															if(current_aic_logP->rank!=index) continue;
															count_gaussians_and_nongaussians(ica_p, distribution, &subgc, &supergc, &gc);
															set_kpvt(ica_p, distribution, kpvt);
															/* Initialize working copy of w */
			
															{
																double *temp=ica_ica_w;
																unsigned int j;
																
																for(j=0; j<ica_p; j++){
																
																	unsigned int i;
																	
																	for(i=0; i<subgc+supergc; i++){
																	
																		if(all_sources==0){
																			if(kpvt[i]==j) *temp++=1.0;
																			else *temp++=0.0;
																		}
																		else{
																			*temp++=null_matrix[kpvt[i]][j];
																		}
																	}
																}
															}
															/* Perform ICA optimization */
															{
																double dx=1.0;
																double eps[2]={1e-4, 1e-8};
																unsigned int maxfun=MAXFUN;
																int stopped=AIR_ucminf(AIR_ica_MLCf15,ica_p*(subgc+supergc),ica_ica_w, &dx, constants, eps, &maxfun, ica_ica_work, ica_p*ica_p*(ica_p*ica_p+13)/2, 1);
																if(stopped < 0){
																	if(subgc+supergc!=0){
#ifndef USING_R
																		printf("%s: %d: ",__FILE__,__LINE__);
																		printf("Failure in ICA routine due to condition %i\n",stopped);
#else
																		REprintf("%s: %d: ",__FILE__,__LINE__);
																		REprintf("Failure in ICA routine due to condition %i\n",stopped);
#endif
																		AIR_free_ica_sources_aic15(constants_ptr);
																		return(AIR_ICA_OPTIMIZATION_ERROR);
																	}
																}
																{
																	// First, incorporate corrections for number of modeled parameters
																	double aic_negloglike=0.0;
																	double negloglike=0;
																	if(working_inf_crit<4){
																		negloglike=AIR_ica_MLCf15(ica_p*ica_p,ica_ica_w,NULL,constants);
																		if(working_inf_crit==1 || working_inf_crit==2){
																			/* AIC and BIC parameter count */
																			aic_negloglike=ica_p*ica_p;
																			aic_negloglike-=(ica_p-gc)*gc;
																			aic_negloglike-=gc*(gc-1)/2;
																			
																			if(!leave_rows_uncentered) aic_negloglike+=ica_p;
																		
																			/* BIC only */
																			if(working_inf_crit==2){
																				aic_negloglike*=log(ica_n)/2.0;
																			}
																		}
																		else if(working_inf_crit==3){
																			if(!fold){
																				/* This is optimized leave-one-out cross validation */
																				ica_n--;
																				{
																					/* Copy ica_ica_w into ica_ica_w2 */
																					double *temp=ica_ica_w;
																					double *temp2=constants[52];
																					
																					unsigned int i;
																					for(i=0; i< ica_p*(subgc+supergc); i++){
																							
																						*temp2++=*temp++;
																					}
																				}
																				{
																					int j;
																					if(samples) j=sample_offset+samples-1;
																					else j=ica_n;
																					for(; j>=0; j--){

																						if(sample_order){
	
																							if(j<sample_offset) break;
																							
																							/* Swap designated column into the last column */
																							{
																								unsigned int i;
																								for(i=0; i<ica_p; i++){
																									
																									double temp=ica_x[sample_order[j]-1][i];
																									ica_x[sample_order[j]-1][i]=ica_x[ica_n][i];
																									ica_x[ica_n][i]=temp;
																								}
																							}
																						}
																						else{
																						
																							/* Swap column j into the last column */
																							unsigned int i;
																							for(i=0; i<ica_p; i++){
																							
																								double temp=ica_x[j][i];
																								ica_x[j][i]=ica_x[ica_n][i];
																								ica_x[ica_n][i]=temp;
																							}
																						}
																						if(!leave_rows_uncentered){
																						
																							/* Re-center based on all but the left out sample
																							 *  and apply same adjustment to the one left out
																							 */
																							unsigned int jj;
																							for(jj=0; jj<=ica_n; jj++){
																								
																								unsigned int i;
																								for(i=0; i<ica_p; i++){
																									
																									double rezero_adj=ica_x[ica_n][i]/ica_n;
																									ica_x[jj][i]+=rezero_adj;
																								}
																							}
																						}
	
																						/* Optimize for all but last column */
																						dx=1.0;
																						maxfun=MAXFUN;
																						stopped=AIR_ucminf(AIR_ica_MLCf15,ica_p*(subgc+supergc),ica_ica_w, &dx, constants, eps, &maxfun, ica_ica_work, ica_p*ica_p*(ica_p*ica_p+13)/2, 1);
																						if(stopped < 0){
																							if(subgc+supergc!=0){
#ifndef USING_R
																								printf("%s: %d: ",__FILE__,__LINE__);
																								printf("Failure in ICA routine due to condition %i\n",stopped);
#else
																								REprintf("%s: %d: ",__FILE__,__LINE__);
																								REprintf("Failure in ICA routine due to condition %i\n",stopped);
#endif
																								ica_n++;
																								AIR_free_ica_sources_aic15(constants_ptr);
																								return(AIR_ICA_OPTIMIZATION_ERROR);
																							}
																						}
																						/* Compute contribution from last column */
																						{
																							double temp=AIR_ica_MLCf15_cv_bias(ica_p*ica_p,ica_ica_w,NULL,constants, 1, &ica_x[ica_n]);
																							aic_negloglike+=temp;
																						}
																						{
																							/* Copy ica_ica_w2 back into ica_ica_w
																							 *  for next iteration or for subsequent
																							 *  computations
																							 */
																							double *temp=ica_ica_w;
																							double *temp2=constants[52];
																							
																							unsigned int i;
																							for(i=0; i< ica_p*(subgc+supergc); i++){
																									
																								*temp++=*temp2++;
																							}
																						}
																						if(sample_order || !leave_rows_uncentered){
																						
																							double **ica_x_copy=constants[12];
																							unsigned int jj;
																							
																							for(jj=0; jj<=ica_n; jj++){
																								
																								unsigned int i;
																								for(i=0; i<ica_p; i++){
																								
																									ica_x[jj][i]=ica_x_copy[jj][i];
																								}
																							}
																						}
																					}
																					if(leave_rows_uncentered && !sample_order){
																						/* Undo the swaps 
																						 * (if leave_rows_uncentered is false or sample_order is not null, original value
																						 * will have just been copied from ica_x_copy
																						 * and would not be recoverable by the
																						 * swaps below
																						 */
																						for(j=0; j<=ica_n; j++){
																							/* Undo the swaps */
																							unsigned int i;
																							for(i=0; i<ica_p; i++){
																								
																								double temp=ica_x[ica_n][i];
																								ica_x[ica_n][i]=ica_x[j][i];
																								ica_x[j][i]=temp;
																							}
																						}
																					}
																				}
																				ica_n++;
																				if(samples) aic_negloglike*=((double)ica_n)/((double)samples);
																			}
																			else{
																				/* This is k-fold cross validation */
																				ica_n-=stride;
																				{
																					/* Copy ica_ica_w into ica_ica_w2 */
																					double *temp=ica_ica_w;
																					double *temp2=constants[52];
																					
																					unsigned int i;
																					for(i=0; i< ica_p*(subgc+supergc); i++){
																							
																						*temp2++=*temp++;
																					}
																				}
																				{
																					/* sample_offset==0, samples=n */
																					int j;
																					for(j=ica_n; j>=0; j-=stride){
																					
																						if(sample_order){
																								
																							/* Swap designated block into the last block */
																							unsigned int s;
																							for(s=0; s<stride; s++){
																							
																								unsigned int i;
																								for(i=0; i<ica_p; i++){
																									
																									double temp=ica_x[sample_order[j+s]-1][i];
																									ica_x[sample_order[j+s]-1][i]=ica_x[ica_n+s][i];
																									ica_x[ica_n+s][i]=temp;
																								}
																							}
																						}
																						else{
																						
																							/* Swap block j into the last block */
																							unsigned int s;
																							for(s=0; s<stride; s++){
																								unsigned int i;
																								for(i=0; i<ica_p; i++){
																								
																									double temp=ica_x[j+s][i];
																									ica_x[j+s][i]=ica_x[ica_n+s][i];
																									ica_x[ica_n+s][i]=temp;
																								}
																							}
																						}
																						if(!leave_rows_uncentered){
																						
																							/* Re-center based on all but the left out block
																							 *  and apply same adjustment to the one left out
																							 */
																							 
																							unsigned int i;
																							for(i=0; i<ica_p; i++){
																							
																								double rezero_adj=0.0;
																								
																								unsigned int jj;
																								for(jj=0; jj<ica_n; jj++){
																									
																									rezero_adj+=ica_x[jj][i];
																								}
																								rezero_adj/=ica_n;
																								for(jj=0; jj<ica_n+stride; jj++){
																								
																									ica_x[jj][i]-=rezero_adj;
																								}
																							}
																						}
	
																						/* Optimize for all but last column */
																						dx=1.0;
																						maxfun=MAXFUN;
																						stopped=AIR_ucminf(AIR_ica_MLCf15,ica_p*(subgc+supergc),ica_ica_w, &dx, constants, eps, &maxfun, ica_ica_work, ica_p*ica_p*(ica_p*ica_p+13)/2, 1);
																						if(stopped < 0){
																							if(subgc+supergc!=0){
#ifndef USING_R
																								printf("%s: %d: ",__FILE__,__LINE__);
																								printf("Failure in ICA routine due to condition %i\n",stopped);
#else
																								REprintf("%s: %d: ",__FILE__,__LINE__);
																								REprintf("Failure in ICA routine due to condition %i\n",stopped);
#endif
																								ica_n++;
																								AIR_free_ica_sources_aic15(constants_ptr);
																								return(AIR_ICA_OPTIMIZATION_ERROR);
																							}
																						}
																						/* Compute contribution from last column */
																						{
																							double temp=AIR_ica_MLCf15_cv_bias(ica_p*ica_p,ica_ica_w,NULL,constants, stride, &ica_x[ica_n]);
																							aic_negloglike+=temp;
																						}
																						{
																							/* Copy ica_ica_w2 back into ica_ica_w
																							 *  for next iteration or for subsequent
																							 *  computations
																							 */
																							double *temp=ica_ica_w;
																							double *temp2=constants[52];
																							
																							unsigned int i;
																							for(i=0; i< ica_p*(subgc+supergc); i++){
																									
																								*temp++=*temp2++;
																							}
																						}
																						if(sample_order || !leave_rows_uncentered || fold){
																						
																							double **ica_x_copy=constants[12];
																							unsigned int jj;
																							
																							for(jj=0; jj<=ica_n; jj++){
																								
																								unsigned int i;
																								for(i=0; i<ica_p; i++){
																								
																									ica_x[jj][i]=ica_x_copy[jj][i];
																								}
																							}
																						}
																					}
																					if(!fold && leave_rows_uncentered && !sample_order){
																						/* Undo the swaps 
																						 * (if fold is non-zero or leave_rows_uncentered is false or samples_order is non-null, original value
																						 * will have just been copied from ica_x_copy
																						 * and would not be recoverable by the
																						 * swaps below
																						 */
																						for(j=0; j<=ica_n; j++){
																							/* Undo the swaps */
																							unsigned int i;
																							for(i=0; i<ica_p; i++){
																								
																								double temp=ica_x[ica_n][i];
																								ica_x[ica_n][i]=ica_x[j][i];
																								ica_x[j][i]=temp;
																							}
																						}
																					}
																				}
																				ica_n+=stride;
																			}
																			aic_negloglike-=negloglike;
																		}
																	}
																	else{
																		/* This uses ica_MLCf15_j2, which has built in bias correction
																		 * for gaussian sources and associated centering terms 
																		 */
																		negloglike=AIR_ica_MLCf15_j2(ica_p*ica_p,ica_ica_w,NULL,constants);
																		/* AIC parameter count for non-gaussians only */
																		aic_negloglike=(ica_p-gc)*ica_p;
																		
																		/* Currently, inf_crit==4 implies !leave_rows_uncentered, but
																		 * test for this for future generality
																		 */
																		if(!leave_rows_uncentered) aic_negloglike+=(ica_p-gc);
																	}
																	// Finally, add computed value
																	aic_negloglike+=negloglike;
																	if(verbose){
#ifndef USING_R
																		printf("cost fxn=%e ", aic_negloglike);
																		printf("subgc supergc gc: %u %u %u ", subgc, supergc, gc);
																		printf("(%u/%u)", index+1, total_permutations);
#else
																		Rprintf("cost fxn=%e ", aic_negloglike);
																		Rprintf("subgc supergc gc: %u %u %u ", subgc, supergc, gc);
																		Rprintf("(%u/%u)", index+1, total_permutations);
#endif
																	}
																	if(aic_negloglike<aic_lowest){
																		if(verbose){
#ifndef USING_R
																			if(pass==passes-1) printf(" best so far\n");
																			else printf("\n");
#else
																			if(pass==passes-1) Rprintf(" best so far\n");
																			else Rprintf("\n");
#endif
																		}
																		aic_lowest=aic_negloglike;
																		{
																			int i;
																			for(i=0;i<ica_p;i++){
																				distribution_aic_best[i]=distribution[i];
																			}
																		}
																	}
																	else{
#ifndef USING_R
																		if(verbose) printf("\n");
#else
																		if(verbose) Rprintf("\n");
#endif
																	}
																	current_aic_logP->value=-aic_negloglike;
																	current_aic_logP->subgc=subgc;
																	current_aic_logP->supergc=supergc;
																	current_aic_logP->gc=gc;
																}
															}
															current_aic_logP++;
														}
													}
													*loglikelihood=-aic_lowest-ica_n*ica_p*log(ica_scale_x);
													{
														unsigned int index;
														for(index=0;index<working_permutations;index++){
															aic_logP[index].value+=aic_lowest;
															aic_logP[index].value=exp(aic_logP[index].value);
														}
													}
													if(passes>1 && pass==0){
														/* Here we find the indices of nondegenerate initialiations */
														/* Sort from highest relative probability to lowest */
														qsort(aic_logP, (size_t)(enough_permutations),sizeof(struct hoods),compare_hoods);
														/* Zero out any item that is degenerate in terms of likelihood and component type */
														if(xval_epsilon<0.0){
															unsigned int i;
															for(i=0; i<enough_permutations; i++){
																char baseval[30];
																char basevaldet[30];
																snprintf(baseval, 29, "%e", aic_logP[i].value);
																snprintf(basevaldet, 29, "%e", aic_logP[i].logdeterminant);
																
																if(aic_logP[i].value==0.0) continue;
																{
																	unsigned int j;
																	for(j=i+1; j<enough_permutations; j++){
																		if(aic_logP[j].subgc==aic_logP[i].subgc && aic_logP[j].supergc==aic_logP[i].supergc && aic_logP[j].gc==aic_logP[i].gc){
																			char newval[30];
																			char newvaldet[30];
																			snprintf(newval, 29, "%e", aic_logP[j].value);
																			snprintf(newvaldet, 29, "%e", aic_logP[j].logdeterminant);
																			if(!strcmp(baseval, newval)){
																				if(!strcmp(basevaldet, newvaldet)){
																					aic_logP[j].value=0.0;
																					aic_logP[j].rank=enough_permutations;
																				}
																			}
																		}
																	}
																}
															}
														}
														else if(xval_epsilon>0.0){
															unsigned int i;
															for(i=0; i<enough_permutations; i++){
																if(aic_logP[i].value==0.0) continue;
																{
																	unsigned int j;
																	for(j=i+1; j<enough_permutations; j++){
																		if(aic_logP[j].subgc==aic_logP[i].subgc && aic_logP[j].supergc==aic_logP[i].supergc && aic_logP[j].gc==aic_logP[i].gc){
																			if(fabs(aic_logP[i].value-aic_logP[j].value)<xval_epsilon){
																				if(fabs(aic_logP[i].logdeterminant-aic_logP[j].logdeterminant)<xval_epsilon){
																					aic_logP[j].value=0.0;
																					aic_logP[j].rank=enough_permutations;
																				}
																			}
																		}
																	}
																}
															}
														}
														/* Sort again based on index used to generate the result */
														qsort(aic_logP, (size_t)(enough_permutations),sizeof(struct hoods),hoods_by_index);
														{
															unsigned int i;
															for(i=0; i<enough_permutations; i++){
																if(aic_logP[i].rank<enough_permutations) working_permutations=i;
															}
														}
														working_permutations++;
													}
												} /* end of passes loop */
												{
													/* Sort from highest relative probability to lowest */
													qsort(aic_logP, (size_t)(working_permutations),sizeof(struct hoods),compare_hoods);
													{
														/* Zero out any item that matches a higher probability item for each type of component */
														unsigned int i;
														for(i=0; i<working_permutations; i++){
														
															if(aic_logP[i].value==0.0) continue;
															{
																unsigned int j;
																for(j=i+1; j<working_permutations; j++){
																	if(aic_logP[j].subgc==aic_logP[i].subgc && aic_logP[j].supergc==aic_logP[i].supergc && aic_logP[j].gc==aic_logP[i].gc){
																		aic_logP[j].value=0.0;
																	}
																}
															}
														}
													}
													/* Sort again from highest relative probability to lowest */
													qsort(aic_logP, (size_t)(working_permutations),sizeof(struct hoods),compare_hoods);
													if(verbose){
#ifndef USING_R
														unsigned int ii;
														printf("\n");
														printf("Relative probabilities of best fitting models for each distinct distributional category:\n");
														for(ii=0; ii<working_permutations; ii++){
															if(aic_logP[ii].value==0.0) break;
															printf("%e subgc=%u supergc=%u gc=%u\n", aic_logP[ii].value, aic_logP[ii].subgc, aic_logP[ii].supergc, aic_logP[ii].gc);
														}
														printf("\n");
#else
														unsigned int ii;
														Rprintf("\n");
														Rprintf("Relative probabilities of best fitting models for each distinct distributional category:\n");
														for(ii=0; ii<working_permutations; ii++){
															if(aic_logP[ii].value==0.0) break;
															Rprintf("%e subgc=%u supergc=%u gc=%u\n", aic_logP[ii].value, aic_logP[ii].subgc, aic_logP[ii].supergc, aic_logP[ii].gc);
														}
														Rprintf("\n");

#endif
													}
												}
												if(model_of_interest!=NULL){
													/* Report on the model of interest */
													*moi_probability=0;
													{
														unsigned int index;
														for(index=0; index<working_permutations; index++){
															if(model_of_interest[0]!=-1){
																if(aic_logP[index].subgc!=(unsigned int)model_of_interest[0]) continue;
															}
															if(model_of_interest[1]!=-1){
																if(aic_logP[index].supergc!=(unsigned int)model_of_interest[1]) continue;
															}
															if(model_of_interest[2]!=-1){
																if(aic_logP[index].gc!=(unsigned int)model_of_interest[2]) continue;
															}
															*moi_probability=aic_logP[index].value;
															break;
														}
													}
												}
												/* Initialize probabilities */
												{
													unsigned int i;
													for(i=0;i<=ica_p;i++){
														probability[i]=0.0;
													}
												}
												{
													/* Find best models with differing numbers of gaussian components */
													unsigned int target_gaussians=min_possible_gaussians;
													
													for(;target_gaussians<=max_possible_gaussians;target_gaussians++){
													
														unsigned int index;
														for(index=0; index<working_permutations; index++){
														
															if(aic_logP[index].gc==target_gaussians){
																probability[target_gaussians]=aic_logP[index].value;
																break;
															}
														}
														if(index==working_permutations) continue; /* no match for this number of gaussians */
														if(aic_logP[index].value==0.0) continue; /* safeguard */
														
														if(outprefix!=NULL && data_saver!=NULL){
														
															if(enough_permutations>1 || inf_crit==3){
																/* Rerun ICA for this item since it is not
																 * guaranteed to be the current working model
																 */
																															
																go_to_index(aic_logP[index].rank, ica_p-nonpermuted_sources, distribution+nonpermuted_sources, permuted_ranges);
																if(distribution[0]==NULL_INDEX){
#ifndef USING_R
																	printf("%s: %d: ",__FILE__,__LINE__);
																	printf("Programming error: failed for find an index previously verified to exist\n");
#else
																	REprintf("%s: %d: ",__FILE__,__LINE__);
																	REprintf("Programming error: failed for find an index previously verified to exist\n");
#endif
																	AIR_free_ica_sources_aic15(constants_ptr);
																	return(AIR_ICA_CODING_ERROR);
																}
																count_gaussians_and_nongaussians(ica_p, distribution, &subgc, &supergc, &gc);
																set_kpvt(ica_p, distribution, kpvt);															
																
																/* Initialize working copy of w */
																{
																	double *temp=ica_ica_w;
																	unsigned int j;
																	
																	for(j=0; j<ica_p; j++){
																	
																		unsigned int i;
																		
																		for(i=0; i<subgc+supergc; i++){
																		
																			if(all_sources==0){
																				if(kpvt[i]==j) *temp++=1.0;
																				else *temp++=0.0;
																			}
																			else{
																				*temp++=null_matrix[kpvt[i]][j];
																			}
																		}
																	}
																}
																{
																	double dx=1.0;
																	double eps[2]={1e-4, 1e-8};
																	unsigned int maxfun=MAXFUN;
																	int stopped=AIR_ucminf(AIR_ica_MLCf15,ica_p*(subgc+supergc),ica_ica_w, &dx, constants, eps, &maxfun, ica_ica_work, ica_p*ica_p*(ica_p*ica_p+13)/2, 1);
																	if(stopped < 0){
																		if(subgc+supergc!=0){
#ifndef USING_R
																			printf("%s: %d: ",__FILE__,__LINE__);
																			printf("Failure in ICA routine due to condition %i\n",stopped);
#else
																			REprintf("%s: %d: ",__FILE__,__LINE__);
																			REprintf("Failure in ICA routine due to condition %i\n",stopped);

#endif
																			AIR_free_ica_sources_aic15(constants_ptr);
																			return(AIR_ICA_OPTIMIZATION_ERROR);
																		}
																	}
																}
																
																/* Force computation of S */
																AIR_ica_MLCf15(ica_p*ica_p, ica_ica_w, NULL, constants);
															}
															{
																/* Globalize variance to that of the original input
																 * matrix
																 */
																unsigned int i;
																
																for(i=0; i<ica_p; i++){
																	variance[i]*=ica_scale_x*ica_scale_x;
																	variance[i]/=total_variance;
																}
															}
															data_saver((AIR_Boolean)FALSE, outprefix, target_gaussians, ica_n, ica_p, ica_s, aic_logP[index].value, 0, subgc, supergc, gc, 0, 0, inf_crit, variance, plot_constants, ow);
														}
													}
												}
							
											}
										}
										if(enough_permutations>1 || inf_crit==3){
										
											/* Must rerun ICA for the globally best model since it is not guaranteed
											 *  to be the current working model
											 */
											{
												/* Reset to the globally best model */
												unsigned int i;
												for(i=0;i<ica_p;i++){
													distribution[i]=distribution_aic_best[i];
												}
												count_gaussians_and_nongaussians(ica_p, distribution, &subgc, &supergc, &gc);
												set_kpvt(ica_p, distribution, kpvt);
											}
											/* Initialize working copy of w */
											{
												double *temp=ica_ica_w;
												unsigned int j;
												
												for(j=0; j<ica_p; j++){
												
													unsigned int i;
													
													for(i=0; i<subgc+supergc; i++){
													
														if(all_sources==0){
															if(kpvt[i]==j) *temp++=1.0;
															else *temp++=0.0;
														}
														else{
															*temp++=null_matrix[kpvt[i]][j];
														}
													}
												}
											}
											
											/* Perform ICA optimization */
											{
												double dx=1.0;
												double eps[2]={1e-4, 1e-8};
												unsigned int maxfun=MAXFUN;
												int stopped=AIR_ucminf(AIR_ica_MLCf15,ica_p*(subgc+supergc),ica_ica_w, &dx, constants, eps, &maxfun, ica_ica_work, ica_p*ica_p*(ica_p*ica_p+13)/2, 1);
												if(stopped < 0){
													if(subgc+supergc!=0){
#ifndef USING_R
														printf("%s: %d: ",__FILE__,__LINE__);
														printf("Failure in ICA routine due to condition %i\n",stopped);
#else
														REprintf("%s: %d: ",__FILE__,__LINE__);
														REprintf("Failure in ICA routine due to condition %i\n",stopped);
#endif
														AIR_free_ica_sources_aic15(constants_ptr);
														return(AIR_ICA_OPTIMIZATION_ERROR);
													}									
												}
											}
											/* Force computation of S */
											AIR_ica_MLCf15(ica_p*ica_p, ica_ica_w, NULL, constants);
											{
												unsigned int i;
												
												for(i=0; i<ica_p; i++){
													variance[i]*=ica_scale_x*ica_scale_x;
													variance[i]/=total_variance;
												}
											}
										}
										else{
											/* The original ICA is the only possible ICA and is still 
											 *  valid (when inf_crit==3, it's not still valid and is
											 *  therefore recalculated).
											 *
											 * variance needs to be globalized to that of the
											 *  original input matrix unless this was already done
											 *  for data_saver routine
											 */
											if(outprefix==NULL || data_saver==NULL){
											
												unsigned int i;

												for(i=0; i<ica_p; i++){
													variance[i]*=ica_scale_x*ica_scale_x;
													variance[i]/=total_variance;
												}
											}
										}
										/* Reorder distribution to correspond to sources, not to mixing matrix */
										{
											unsigned int i;
											for(i=0; i<subgc; i++){
												distribution[i]=0;
											}
											for(; i<subgc+supergc; i++){
												distribution[i]=1;
											}
											for(; i<ica_p; i++){
												distribution[i]=2;
											}
										}
									}
								}
								*ica_s_ptr=ica_s;				/* Don't free ica_s after this */
								*distribution_ptr=distribution;	/* Don't free distribution after this */
								*variance_ptr=variance;			/* Don't free variance after this */
								*probability_ptr=probability;	/* Don't free probability after this */
							}
						}	
					}
				}
			}
		}
	}
	return 0;
 }

