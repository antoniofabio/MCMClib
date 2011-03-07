/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
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

/** \brief alloc a new AT7 sampler
@param r RNG
@param logdistr target log-distrib. fun.
@param logdistr_data log-distrib. fun. extra data
@param x current chain value
@param t0 burn in length
@param beta_hat starting mixture weights estimates
@param mu_hat starting mixture means estimates
@param Sigma_hat starting mixture variance estimates
*/
mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p logdistr, void* logdistr_data,
			       gsl_vector* x, size_t t0,
			       const gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat);

/**\brief set scaling factors of an at7 sampler */
void mcmclib_at7_set_sf(mcmclib_amh* p, const gsl_vector* sf);
/**\brief set all scaling factors of an at7 sampler to the same value */
void mcmclib_at7_set_sf_all(mcmclib_amh* p, const double all);

/**@}*/
/**@}*/
#endif
