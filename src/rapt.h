/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __RAPT_H__
#define __RAPT_H__

#include "amh.h"
#include "rapt_q.h"

/**\addtogroup adaptive
@{*/
/**\defgroup RAPT rapt
\brief Regional AdaPTive
@{*/

/** alloc a new RAPT sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x current chain value
@param t0 burn-in length before starting adaptation
@param sigma_whole global proposal covariance matrix
@param K number of regions
@param sigma_local array of local proposal covariance matrices
@param which_region boundary computing function
@param which_region_data ptr to extra 'which_region' data
@param region_data_free 'which_region' func data destructor
*/
mcmclib_amh* mcmclib_rapt_alloc(gsl_rng* r,
				distrfun_p logdistr, void* logdistr_data,
				gsl_vector* x,
				size_t t0,
				const gsl_matrix* sigma_whole,
				size_t K,
				gsl_matrix** sigma_local,
				region_fun_t which_region,
				void* which_region_data,
				free_fun_t region_data_free);

/** customly set additive variance correction factor */
void mcmclib_rapt_set_correction_factor(mcmclib_amh* p, const double eps);

/** customly set local and global scaling factors */
void mcmclib_rapt_set_scaling_factors(mcmclib_amh* p, const double local, const double global);

/**@}*/
/**@}*/
#endif
