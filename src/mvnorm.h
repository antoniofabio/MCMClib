/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MVNORM_H__
#define __MCMCLIB_MVNORM_H__

#include "common.h"
/**\addtogroup distributions
 @{*/

/**\addtogroup multivariate
 @{*/

/** Multivariate gaussian variate
@param r RNG state
@param sigma variance/covariance matrix
@param out result vector
*/
void mcmclib_mvnorm(const gsl_rng* r,
		    const gsl_matrix* sigma,
		    gsl_vector* out);

/** Multivariate gassian variate, with known cholesky decomposition.
\internal
@param r RNG state
@param sigma_chol cholesky decomposition of var/covariance matrix
@param out result vector
*/
void mcmclib_mvnorm_chol(const gsl_rng* r,
	const gsl_matrix* sigma_chol,
	gsl_vector* out);

/** Multivariate gassian variate, with known cholesky dec. of the precision matrix.
@param r RNG state
@param Psi cholesky decomposition of the precision matrix
@param out result vector
*/
void mcmclib_mvnorm_precision(const gsl_rng* r,
			      const gsl_matrix* Psi,
			      gsl_vector* out);

/**\brief Multivariate Gaussian distribution*/
typedef struct {
  gsl_vector* mean; /**< distribution mean*/
  double* vcov; /**< distribution variance/covariance matrix*/
  gsl_matrix* rooti;
  gsl_vector* x_mu;
  gsl_vector* mahal;
  double determinant;
} mcmclib_mvnorm_lpdf;

/** Alloc extra data for a multivariate gaussian distribution
@param mean mean
@param vcov variance/covariance matrix (raw data)
*/
mcmclib_mvnorm_lpdf* mcmclib_mvnorm_lpdf_alloc(gsl_vector* mean, double* vcov);

/** Free extra data for a multivariate gaussian distribution
@param p pointer to distrib extra data
*/
void mcmclib_mvnorm_lpdf_free(mcmclib_mvnorm_lpdf* p);

/** Multivariate gassian log-distribution
@param in_p extra data, allocated via \ref mcmclib_mvnorm_lpdf_alloc
@return log-pdf
*/
double mcmclib_mvnorm_lpdf_compute(void* in_p, gsl_vector* x);

/**update cholesky decomposition info
   @return !0 if non pos.def. covariance matrix
\internal*/
int mcmclib_mvnorm_lpdf_chol(mcmclib_mvnorm_lpdf* p);

/**compute log-distrib without recomputing cholesky decomposition
\internal*/
double mcmclib_mvnorm_lpdf_compute_nochol(mcmclib_mvnorm_lpdf* p, gsl_vector* x);

/**update inverse covariance matrix info
\internal*/
void mcmclib_mvnorm_lpdf_inverse(mcmclib_mvnorm_lpdf* p);
/**compute log-distrib by exploiting previously computed inverse
\internal*/
double mcmclib_mvnorm_lpdf_compute_noinv(mcmclib_mvnorm_lpdf* p, gsl_vector* x);

/** Multivariate normal distribution from precision matrix
    @param mu mean
    @param iSigma precision matrix
    @param x
    @param ldet log-determinant of var/cov matrix
    @param work1 workspace vector of size equal to that of 'x'
    @param work2 same as work1
    @return log-pdf
 */
double mcmclib_mvnorm_lpdf_noinv(gsl_vector* mu, gsl_matrix* iSigma, gsl_vector* x,
				 double ldet, gsl_vector* work1, gsl_vector* work2);

/** multivariate zero-mean normal log-density based on precision matrix
    
    The lower triangular part of 'Psi' is referenced. The upper triangular part
    is ignored. */
double mcmclib_mvnormzp_lpdf(const gsl_matrix* Psi, const gsl_vector* y);

/**@}*/
/**@}*/
#endif
