/* Copyright 2010-2014 Roger P. Woods, M.D. */
/* Modified: 10/19/2014 */

#include "AIR.h"

static void ica_pca_free(void *ptr){
#ifndef USING_R
	free(ptr);
#endif
}
static void *ica_pca_malloc(size_t n, int size){
#ifndef USING_R
	return(malloc(n*size));
#else
	return(R_alloc(n, size));
#endif
}

void ica_pca_R(int *error, int *seed, int *offset_random, 
int *sample_order, int *sample_order_length, int *sample_offset, int *samples,
int *inf_crit, int *fold, double *xval_epsilon, 
int *ranges, double *em_array, int *mrows, int *ncols, int *components, 
double *ica_s_array, double *loglikelihood, int *desired_initialization,
int *distribution_R, double *variance_R, 
double *probability_R, int *hinted_subgaussian_sources, 
int *hinted_supergaussian_sources, int *hinted_unspecified_sources,
double *source_array, int *leave_rows_uncentered
){

	void **constants=NULL;
	double **ica_s=NULL;
	unsigned int *distribution=NULL;
	double *variance=NULL;
	double *probability=NULL;
	double **source_matrix=NULL;
		
	if(*sample_order_length==0){
		sample_order=NULL;
	}
			
	unsigned int source_count=*hinted_subgaussian_sources+*hinted_supergaussian_sources+*hinted_unspecified_sources;
	double **em=(double **)ica_pca_malloc((size_t)*ncols, sizeof(double *));
	if(!em){
		*error=1;
		return;
	}
	{
		unsigned int j;
		for(j=0; j<*ncols; j++){
			em[j]=em_array+j**mrows;
		}
	}
	if(source_count>0){
		source_matrix=(double **)ica_pca_malloc((size_t)source_count, sizeof(double *));
		if(!source_matrix){
			ica_pca_free(em);
			*error=1;
			return;
		}
		{
			unsigned int j;
			for(j=0; j<source_count; j++){
				source_matrix[j]=source_array+j**ncols;
			}
		}
	}
	{
		AIR_Error errcode=AIR_ica_sources_aic15(*seed, *offset_random, (unsigned int *)sample_order, *sample_offset, *samples, &constants, *inf_crit, *fold, *xval_epsilon, (unsigned int *)ranges, em, *mrows, *ncols, *hinted_subgaussian_sources, *hinted_supergaussian_sources, *hinted_unspecified_sources, source_matrix, (unsigned int *)components, &ica_s, *desired_initialization, &distribution, &variance, &probability, loglikelihood, NULL, NULL, NULL, NULL, NULL, (AIR_Boolean)*leave_rows_uncentered, (AIR_Boolean)FALSE, (AIR_Boolean)FALSE);
		if(errcode){
			ica_pca_free(em);
			if(source_matrix!=NULL) ica_pca_free(source_matrix);
			AIR_free_ica_sources_aic15(&constants);
			*error=1;
			return;
		}
	}
	{
		unsigned int j;
		for(j=0; j<*ncols; j++){
			unsigned int i;
			for(i=0; i<*components; i++){
				*ica_s_array++=ica_s[j][i];
			}
		}
	}
	
	{
		unsigned int i;
		for(i=0; i<*components; i++){
			*distribution_R++=distribution[i];
			*variance_R++=variance[i];
			*probability_R++=probability[i];
		}
		*probability_R++=probability[i];
		
	}
	
	ica_pca_free(em);
	if(source_matrix!=NULL) ica_pca_free(source_matrix);
	AIR_free_ica_sources_aic15(&constants);
	*error=0;
}

