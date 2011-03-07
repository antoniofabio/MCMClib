/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __GAUSS_SCALAR_AM_H__
#define __GAUSS_SCALAR_AM_H__

/**\addtogroup adaptive
@{
\defgroup GAUSS_SCALAR_AM Scalar Gaussian random walk Adaptive Metropolis
@{*/

#include "amh.h"

/** alloc (and init) extra AM data
@param r RNG state
@param logdistr pointer to a log-distribution function
@param logdistr_data extra data to be passed to the log-distribution function
@param x starting value
@param scaling proposal variance scaling factor
@param N0 burn-in length before starting adaptation
*/
mcmclib_amh* mcmclib_gauss_scalar_am_alloc(gsl_rng* r, distrfun_p logdistr, void* logdistr_data,
					   gsl_vector* x, const double scaling, const size_t N0);

/**@}*/
/**@}*/
#endif
