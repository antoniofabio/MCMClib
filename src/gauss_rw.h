/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include "mh.h"

/**\addtogroup metropolis_samplers
@{
\defgroup gauss_rw Component-wise Gaussian Random Walk
@{*/

void mcmclib_gauss_rw_sample(mcmclib_mh_q* q, gsl_vector* x);

/** alloc (and init) Gaussian RW object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param start_x starting value
@param data extra data to be passed to the distribution function
@param step_size gaussian proposal width (s.d.)
*/
mcmclib_mh* mcmclib_gauss_rw_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* data,
				   gsl_vector* start_x, double step_size);

/**@}*/
/**@}*/
#endif
