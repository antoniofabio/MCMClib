/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_LM_H__
#define __MCMCLIB_LM_H__

#include "mvnorm.h"
/**\addtogroup Models
 @{*/

/**\addtogroup LinearModel
 @{*/

/**\brief Linear regression model with conjugte priors */
typedef struct {
} mcmclib_lm;

/** Alloc extra data for an lm model
    @param X design matrix
 */
mcmclib_lm* mcmclib_lm_alloc(const gsl_matrix* X);

/** Free extra data for an lm object */
void mcmclib_lm_free(mcmclib_lm* p);


/** Linear regression with Gaussian errors beta draw
    (Normal  prior).

    Ported from the MCMCpack R package 'NormNormregress_beta_draw' C++ sources

    regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)

    @param rng RNG
    @param XX X'X product
    @param XY X'Y product
    @param b0 prior mean for beta
    @param B0 prior var/cov matrix for beta
    @param sigma2 known errors variance
    @param out sampled value goes here
    @return GSL_SUCCESS if all goes right
*/
int mcmclib_lm_beta_draw(gsl_rng* rng, const gsl_matrix* XX, const gsl_vector* XY,
			 const gsl_vector* b0, const gsl_matrix* B0,
			 double sigma2, gsl_vector* out);

/** Linear regression with Gaussian errors sigma2 draw
    (inverse-Gamma  prior).

    Ported from the MCMCpack R package C++ sources.

    regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)

    @param c0 c0/2 is the prior shape parameter for sigma2
    @param d0 d0/2 is the prior scale parameter for sigma2
    @return GSL_SUCCESS if all goes right
*/
int mcmclib_lm_sigma2_draw (gsl_rng* rng, const gsl_matrix* X, const gsl_vector* Y,
			    const gsl_vector* beta, double c0, double d0,
			    gsl_vector* out);

/**@}*/
/**@}*/
#endif
