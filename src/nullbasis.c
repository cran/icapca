/* Copyright 2004-2010 Roger P. Woods, M.D. */
/* Modified 11/19/2010 */

/*
 *	Given an input matrix INPUT[COLUMNS][ROWS] with ROWS>COLUMNS, this will return, in the
 *  last (ROWS-COLUMNS) columns of OUTPUT, a basis for the null space not spanned by the vectors
 *	that make up the columns of INPUT. QR decomposition is used for computation.
 *
 *	If JPVT is not NULL, it is assumed to be initialized as required for dqrdc. Pivoting will
 * 	only be performed in dqrdc if JPVT is not null
 *
 *	If WORK and/or QRAUX are NULL, they will be appropriately allocated and freed internally
 *	Note that WORK will not be used if JPVT is null
 *
 *
 *	On Entry:
 *		input[columns][rows]
 *		rows>columns
 *		output[rows][rows]
 *		jpvt[columns] or NULL
 *		work[columns] or NULL (unused if jvpt is NULL)
 *		qraux[columns] or NULL
 *
 *	On Exit:
 *		input will have been modified to contain the output of dqrdc
 *		jpvt, if not null, will have been modified to reflect pivoting as set in dqrdc
 *		output[rows][rows] will have the null space basis in its final (rows-columns) columns
 *							The initial columns of output will contain an orthonormal basis for
 *							the subspace spanned by the columns of the original input
 *
 *	See dqrdc for more information regarding pivoting when input is rank deficient
 */

#include "AIR.h"

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
	return(malloc(n*size));
#else
	return(mxMalloc(n*size));
#endif
#else
	return(R_alloc(n, size));
#endif
}
static double *ica_pca_matrix1(const size_t a)
{
	double *high=(double *)ica_pca_malloc(a, sizeof(double));
	if(!high) return NULL;
	return high;
}

static void free_function(AIR_Boolean free_work, double *work, AIR_Boolean free_qraux, double *qraux)
{
	if(free_work) ica_pca_free(work);
	if(free_qraux) ica_pca_free(qraux);
}

AIR_Error AIR_nullbasis(double **input, const unsigned int rows, const unsigned int columns, double **output, int *jpvt, double *work, double *qraux)
{
	AIR_Boolean free_work=FALSE;
	AIR_Boolean free_qraux=FALSE;
	
	if(jpvt!=NULL){
		if(work==NULL){
			work=ica_pca_matrix1(columns);
			if(!work){
				free_function(free_work, work, free_qraux, qraux);
				return AIR_MEMORY_ALLOCATION_ERROR;
			}
			else{
				free_work=TRUE;
			}
		}
	}
	if(qraux==NULL){
		qraux=ica_pca_matrix1(columns);
		if(!qraux){
			free_function(free_work, work, free_qraux, qraux);
			return AIR_MEMORY_ALLOCATION_ERROR;
		}
		else{
			free_qraux=TRUE;
		}
	}
	AIR_dqrdc(input, rows, columns, jpvt, work, jpvt!=NULL, qraux);
	{
		unsigned int j;
		for(j=0;j<rows;j++){
			unsigned int i;
			for(i=0;i<rows;i++){
				if(i==j) output[j][i]=1.0;
				else output[j][i]=0.0;
			}
			{
				int info;	// No need to check info afterwards when job==10000
				AIR_dqrsl(input, rows, columns, qraux, output[j], 10000, output[j], NULL, NULL, NULL, NULL, &info);
			}
		}
	}
	free_function(free_work, work, free_qraux, qraux);
	return 0;
}
