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

/** GRW proposal log-density (fake) \internal*/
double mcmclib_gauss_mrw_qd(mcmclib_mh_q* ignore, gsl_vector* x, gsl_vector* y);
/** GRW proposal sampler */
void mcmclib_gauss_mrw_sample(mcmclib_mh_q* q, gsl_vector* x);

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


/**@}*/
/**@}*/
#endif
