/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __GAUSS_AM_H__
#define __GAUSS_AM_H__

/**\addtogroup adaptive
@{
\defgroup GAUSS_AM Gaussian random walk Adaptive Metropolis
@{*/

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "amh.h"
#include "gauss_mrw.h"

/**\brief Gaussian AM cumulated sufficient statistics and support data*/
typedef struct {
  gsl_vector* sum_x; /**< cumulated sum of xs*/
  gsl_matrix* sum_xx; /**< cumulated sum of xxs*/
  gsl_matrix* Sigma_eps; /**< pos. definiteness cov. correction additive constant*/
  gsl_matrix* Sigma_zero; /**< starting proposal covariance matrix*/
  size_t t0; /**< burn in before starting adaptation*/
  double sf; /**< scaling factor*/
} mcmclib_gauss_am_suff;

/** free extra AM data*/
void mcmclib_gauss_am_suff_free(void* p);

/** alloc (and init) extra AM data
@param r RNG state
@param logdistr pointer to a log-distribution function
@param logdistrib_data extra data to be passed to the log-distribution function
@param start_x starting value
@param sigma_zero starting proposal covariance matrix
@param t0 burn-in length before starting adaptation
*/
mcmclib_amh* mcmclib_gauss_am_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* logdistr_data,
				    gsl_vector* start_x,
				    const gsl_matrix* sigma_zero, size_t t0);

/** set scaling factor for the gaussian AM sampler 'p'*/
void mcmclib_gauss_am_set_sf(mcmclib_amh* p, double sf);

/** AM gamma update function \internal */
void mcmclib_gauss_am_update_gamma(mcmclib_amh* p);

/**@}*/
/**@}*/
#endif
