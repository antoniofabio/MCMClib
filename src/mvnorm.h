#ifndef __MCMCLIB_MVNORM_H__
#define __MCMCLIB_MVNORM_H__

#include "common.h"

/** multivariate gaussian variate
@param r RNG state
@param sigma variance/covariance matrix
@param out result vector
*/
void mcmclib_mvnorm(const gsl_rng* r,
	const gsl_matrix* sigma,
	gsl_vector* out);

/** multivariate gassian variate, with known cholesky decomposition
@param r RNG state
@param sigma_chol cholesky decomposition of var/covariance matrix
@param out result vector
*/
void mcmclib_mvnorm_chol(const gsl_rng* r,
	const gsl_matrix* sigma_chol,
	gsl_vector* out);

#endif
