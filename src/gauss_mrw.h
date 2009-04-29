/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __GAUSS_MRW_H__
#define __GAUSS_MRW_H__

/**\addtogroup metropolis_samplers
@{
\defgroup GAUSS_MRW Multivariate Gaussian Random Walk
@{*/
#include "mh.h"

typedef struct {
  gsl_rng* r;
  gsl_matrix* Sigma;
} mcmclib_gauss_mrw_gamma;

mcmclib_gauss_mrw_gamma* mcmclib_gauss_mrw_gamma_alloc(gsl_rng* r, const gsl_matrix* Sigma);
void mcmclib_gauss_mrw_gamma_free(mcmclib_gauss_mrw_gamma* p);
void mcmclib_gauss_mrw_sample(void* in_p, gsl_vector* x);
double mcmclib_gauss_mrw_qd(void* ignore, gsl_vector* x, gsl_vector* y);

mcmclib_mh_q* mcmclib_gauss_mrw_q_alloc(gsl_rng* r, const gsl_matrix* Sigma);
void mcmclib_gauss_mrw_q_free(mcmclib_mh_q* p);

/** alloc (and init) a Gaussian RW object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param data extra data to be passed to the distribution function
@param start_x starting value
@param sigma_prop gaussian proposal covariance matrix
*/
mcmclib_mh* mcmclib_gauss_mrw_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* data,
				    gsl_vector* start_x, const gsl_matrix* sigma_prop);

/** free extra Gaussian RW data*/
void mcmclib_gauss_mrw_free(mcmclib_mh* p);

/** GRW proposal log-density (fake) \internal*/
double mcmclib_gauss_mrw_qd(void* ignore, gsl_vector* x, gsl_vector* y);

/** GRW proposal sampler \internal
@param in_p ptr to a variance-covariance matrix
@param x object in which to put result
*/
void mcmclib_gauss_mrw_sample(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
