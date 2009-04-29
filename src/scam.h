/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __SCAM_H__
#define __SCAM_H__

/**\addtogroup adaptive
@{
\defgroup SCAM Single Component Adaptive Metropolis
@{*/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#include "amh.h"

/** pointer to a vector of distribution functions */
typedef double (*distrfun_i_p) (void* data, int index, double x);

/**\brief SCAM support data*/
typedef struct {
	gsl_vector* x_full; /**< full state vector*/
	gsl_vector* xi; /**< univariate state view*/
	mcmclib_amh** x_smp; /**< array of univariate adaptive samplers*/
	distrfun_i_p logdistr; /**< target log-distribution*/
	void* logdistr_data; /**< target log-distribution extra data*/
	int curr_index; /**< utility field to pass extra argument to user logdistrib function*/
} mcmclib_scam;

/** alloc (and init) SCAM data
@param r RNG state
@param logdistr pointer to a log-distribution function
@param logdistrib_data extra data to be passed to the log-distribution function
@param x full state vector
@param sigma_zero starting proposal variance
@param t0 burn-in length before starting adaptation
*/
mcmclib_scam* mcmclib_scam_alloc(gsl_rng* r,
																 distrfun_i_p logdistr, void* logdistr_data,
																 gsl_vector* x,
																 double sigma_zero, int t0);
/** free SCAM data*/
void mcmclib_scam_free(mcmclib_scam* p);
/** update all state components, one after another */
void mcmclib_scam_update(mcmclib_scam* p);

/** SCAM target logdistrib wrapper function \internal*/
double mcmclib_scam_logdistr(void* data, gsl_vector* x);

/**@}*/
/**@}*/
#endif