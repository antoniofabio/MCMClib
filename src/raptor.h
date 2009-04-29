/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __RAPTOR_H__
#define __RAPTOR_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "mixem_online.h"
#include "mixnorm.h"
#include "rapt.h"

/**\addtogroup adaptive
@{*/
/**\defgroup RAPTOR raptor
\brief RAPT based on On-Line EM fitting of a Gaussian mixture
@{*/

/** \brief RAPTOR sampler gamma values */
typedef struct {
  gsl_vector* beta_hat; /**< current mixture weights estimates*/
  gsl_vector** mu_hat; /**< current mixture means estimates*/
  gsl_matrix** Sigma_hat; /**< current mixture variances estimates*/

  mcmclib_mvnorm_lpdf** pik_hat; /**< single mixture components densities*/
  mcmclib_mixnorm_lpdf* pi_hat; /**< mixture density*/
} mcmclib_raptor_gamma;

/** \brief RAPTOR sufficient data */
typedef struct {
  mcmclib_mixem_online* em; /**< online-EM mixture fitter*/

  gsl_matrix* Sigma_eps;
  double scaling_factor_local;
  double scaling_factor_global;
} mcmclib_raptor_suff;

/** alloc a new RAPTOR sampler suff. stats. object
@param t0 burn-in length before starting adaptation
@returns a new raptor_suff object
*/
mcmclib_raptor_suff* mcmclib_raptor_suff_alloc(mcmclib_raptor_gamma* g, int t0);
/** free raptor_suff data*/
void mcmclib_raptor_suff_free(mcmclib_raptor_suff* p);
/** Update suff. stats. of a RAPTOR chain*/
int mcmclib_raptor_suff_update(mcmclib_raptor_suff* p);

/** \brief alloc a new RAPTOR sampler
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
mcmclib_amh* mcmclib_raptor_alloc(gsl_rng* r,
				  distrfun_p logdistr, void* logdistr_data,
				  gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
				  gsl_vector* beta_hat,
				  gsl_vector** mu_hat,
				  gsl_matrix** Sigma_hat);
/**\brief free a previously allocated RAPTOR sampler*/
void mcmclib_raptor_free(mcmclib_amh* p);

/** update local and global RAPT proposals covariance matrices

basing on current mixture parameters estimates*/
void mcmclib_raptor_update_proposals(mcmclib_raptor_suff* p);

/**@}*/
/**@}*/
#endif
