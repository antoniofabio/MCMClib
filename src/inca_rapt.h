/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __INCA_RAPT_H__
#define __INCA_RAPT_H__

#include "inca.h"
#include "rapt.h"

/**\addtogroup adaptive
@{*/
/**\defgroup INCA_RAPT inca_rapt
\brief INCA Regional AdaPTive
@{*/

/** alloc a new INCA RAPT sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x array of current chain values
@param M number of parallel chains
@param t0 burn-in length before starting adaptation
@param sigma_whole global proposal covariance matrix
@param K number of regions
@param sigma_local array of local proposal covariance matrices
@param which_region boundary computing function
@param which_region_data ptr to extra 'which_region' data
*/
mcmclib_inca* mcmclib_inca_rapt_alloc(gsl_rng* r,
				      distrfun_p logdistr, void* logdistr_data,
				      gsl_vector** x, int M, int t0,
				      const gsl_matrix* sigma_whole,
				      int K, gsl_matrix** sigma_local,
				      region_fun_t which_region,
				      void* which_region_data);

/** free  data*/
void mcmclib_inca_rapt_free(mcmclib_inca* p);

/**@}*/
/**@}*/
#endif
