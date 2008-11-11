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

/** multivariate gaussian distribution parameters
*/
typedef struct {
  gsl_vector* mean;
  double* vcov;
  gsl_matrix* rooti;
  gsl_vector* x_mu;
  gsl_vector* mahal;
  double determinant;
} mcmclib_mvnorm_lpdf;

/** alloc extra data for a multivariate gaussian distribution
@param mean mean
@param vcov variance/covariance matrix (raw data)
*/
mcmclib_mvnorm_lpdf* mcmclib_mvnorm_lpdf_alloc(gsl_vector* mean, double* vcov);

/** free extra data for a multivariate gaussian distribution
@param p pointer to distrib extra data
*/
void mcmclib_mvnorm_lpdf_free(mcmclib_mvnorm_lpdf* p);

/** multivariate gassian log-distribution
@param in_p extra data, allocated via 'mcmclib_mvnorm_lpdf_alloc'
@return log-pdf
*/
double mcmclib_mvnorm_lpdf_compute(void* in_p, gsl_vector* x);

/**update cholesky decomposition info*/
void mcmclib_mvnorm_lpdf_chol(mcmclib_mvnorm_lpdf* p);

/**compute log-distrib without recomputing cholesky decomposition*/
double mcmclib_mvnorm_lpdf_compute_nochol(mcmclib_mvnorm_lpdf* p, gsl_vector* x);

/**update inverse covariance matrix info*/
void mcmclib_mvnorm_lpdf_inverse(mcmclib_mvnorm_lpdf* p);
/**compute log-distrib by exploiting previously computed inverse*/
double mcmclib_mvnorm_lpdf_compute_noinv(mcmclib_mvnorm_lpdf* p, gsl_vector* x);

#endif
