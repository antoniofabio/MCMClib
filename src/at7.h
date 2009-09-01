/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __AT7_H__
#define __AT7_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "mixem_online.h"
#include "mixnorm.h"
#include "raptor.h"

/**\addtogroup adaptive
@{*/
/**\defgroup AT7 AT7
\brief Andrieu and Thoms (2008) Alg. 7
@{*/

/** \brief AT7 gamma values */
typedef struct {
  gsl_vector* beta_hat; /**< current mixture weights estimates*/
  gsl_vector** mu_hat; /**< current mixture means estimates*/
  gsl_matrix** Sigma_hat; /**< current mixture variances estimates*/

  mcmclib_mvnorm_lpdf** pik_hat; /**< single mixture components densities*/
  mcmclib_mixnorm_lpdf* pi_hat; /**< mixture density*/

  gsl_matrix* Sigma_eps; /**< positive-definiteness correction factor*/
  gsl_vector* scaling_factors; /**< region-specific scaling factors*/
} mcmclib_at7_gamma;

/** alloc a new at7_gamma object. Input arguments are copied @internal */
mcmclib_at7_gamma* mcmclib_at7_gamma_alloc(gsl_vector* beta_hat,
					   gsl_vector** mu_hat,
					   gsl_matrix** Sigma_hat);
/** frees an at7_gamma object @internal*/
void mcmclib_at7_gamma_free(mcmclib_at7_gamma* p);

/** \brief AT7 sufficient data */
typedef struct {
  mcmclib_mixem_online* em; /**< online-EM mixture fitter*/
} mcmclib_at7_suff;

/** alloc a new AT7 sampler suff stats object
@param t0 burn-in length before starting adaptation
@returns a new AT7_suff object
*/
mcmclib_at7_suff* mcmclib_at7_suff_alloc(mcmclib_at7_gamma* g, int t0,
					 mcmclib_rapt_gamma* rg);
/** free raptor_suff data*/
void mcmclib_at7_suff_free(mcmclib_at7_suff* p);
/** Update suff stats of an AT7 chain*/
int mcmclib_at7_suff_update(mcmclib_raptor_suff* p);

/** \brief alloc a new AT7 sampler
@param r RNG
@param logdistr target log-distrib. fun.
@param logditsr_data log-distrib. fun. extra data
@param x current chain value
@param t0 burn in length
@param Sigma_zero starting variance-covariance matrix
@param beta_hat starting mixture weights estimates
@param mu_hat starting mixture means estimates
@param Sigma_hat starting mixture variance estimates
*/
mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p logdistr, void* logdistr_data,
			       gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
			       gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat);
/**\brief free a previously allocated AT7 sampler*/
void mcmclib_at7_free(mcmclib_amh* p);

/**\brief set scaling factors */
void mcmclib_at7_set_sf(mcmclib_amh* p, const gsl_vector* sf);
/**\brief set all scaling factors to the same value */
void mcmclib_at7_set_sf_all(mcmclib_amh* p, double all);

/**@}*/
/**@}*/
#endif
