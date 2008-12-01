#ifndef __VECTOR_STATS_H__
#define __VECTOR_STATS_H__
#include <gsl/gsl_statistics_double.h>
#include "common.h"

/**\addtogroup descriptive
@{*/

/** column means */
void mcmclib_matrix_colmeans(gsl_matrix* m, gsl_vector* out);
/** row means */
void mcmclib_matrix_rowmeans(gsl_matrix* m, gsl_vector* out);

/** get variance/covariance matrix out of the 'vertical' matrix 'm' */
void mcmclib_matrix_covariance(gsl_matrix* m, gsl_matrix* out);

/** update covariance value 'recursively'
@param cov current covariance matrix
@param mean current mean
@param n current sample size
@param x new value
@return nothing. 'cov', 'mean' and 'n' values will be updated as a side-effect
*/
void mcmclib_covariance_update(gsl_matrix* cov, gsl_vector* mean, int* n, gsl_vector* x);

/**Pooled weighted variance
@param beta first component's weight
@param means array of 2 means
@param variances array of two variances
@param output result
*/
void mcmclib_pooled_variance(double beta, gsl_vector** means,
			     gsl_matrix** variances, gsl_matrix* V);

/**@}*/

#endif
