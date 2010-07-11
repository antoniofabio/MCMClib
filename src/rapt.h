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

/**\brief RAPT sufficient statistics*/
typedef struct {
  size_t t0; /**< burn-in length*/
  gsl_vector** means; /**< array of regions means*/
  gsl_matrix** variances; /**< array of regions variances*/
  gsl_vector* global_mean;
  gsl_matrix* global_variance;
  gsl_vector* n; /**< number of visits in each region*/

  /*internal data*/
  gsl_matrix* Sigma_eps; /**< additive perturbation factor for variances updating*/
  gsl_vector* workspace; /**< utility workspace memory*/
  double scaling_factor_local; /**< local proposal variance scaling factor*/
  double scaling_factor_global; /**< global proposal variance scaling factor*/
} mcmclib_rapt_suff;

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
*/
mcmclib_amh* mcmclib_rapt_alloc(gsl_rng* r,
				distrfun_p logdistr, void* logdistr_data,
				gsl_vector* x,
				size_t t0,
				const gsl_matrix* sigma_whole,
				size_t K,
				gsl_matrix** sigma_local,
				region_fun_t which_region,
				void* which_region_data);

/** free  data*/
void mcmclib_rapt_free(mcmclib_amh* p);

/** update local and global proposals covariance matrices
basing on current (region-specific) sample variances*/
void mcmclib_rapt_update_proposals(mcmclib_amh* p);

/** update local and global proposals covariance matrices
using custom estimated local and global variances
 */
void mcmclib_rapt_update_proposals_custom(mcmclib_amh* p,
					  gsl_matrix** variances,
					  gsl_matrix* global_variance);

/** customly set additive variance correction factor */
void mcmclib_rapt_suff_set_correction_factor(mcmclib_rapt_suff* p, double eps);

/** customly set local and global scaling factors */
void mcmclib_rapt_suff_set_scaling_factors(mcmclib_rapt_suff* p, double local, double global);

/**@}*/
/**@}*/
#endif
