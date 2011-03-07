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

/** \brief alloc a new RAPTOR sampler
@param r RNG
@param logdistr target log-distrib. fun.
@param logdistr_data log-distrib. fun. extra data
@param x current chain value
@param t0 burn in length
@param Sigma_zero starting variance-covariance matrix
@param beta_hat starting mixture weights estimates
@param mu_hat starting mixture means estimates
@param Sigma_hat starting mixture variance estimates
*/
mcmclib_amh* mcmclib_raptor_alloc(gsl_rng* r,
				  distrfun_p logdistr, void* logdistr_data,
				  gsl_vector* x, size_t t0, gsl_matrix* Sigma_zero,
				  const gsl_vector* beta_hat,
				  gsl_vector** mu_hat,
				  gsl_matrix** Sigma_hat);

/**\brief set scaling factor */
void mcmclib_raptor_set_sf(mcmclib_amh* p, double sf);
/**\brief set global scaling factor */
void mcmclib_raptor_set_sf_global(mcmclib_amh* p, double sf);
/**\brief set local scaling factor */
void mcmclib_raptor_set_sf_local(mcmclib_amh* p, double sf);
/** customly set global proposal weight (same for all regions)*/
void mcmclib_raptor_set_alpha(mcmclib_amh* p, double alpha);

/**@}*/
/**@}*/
#endif
