/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MCAR_MODEL_H__
#define __MCMCLIB_MCAR_MODEL_H__
/**\addtogroup Models
 @{*/
/**\addtogroup MCAR
 @{*/

#include "mcar_tilde.h"
#include "lpdf_iwishart.h"

/**\brief MCAR(\tilde B, \Gamma) Model */
typedef struct {
  mcmclib_mcar_tilde_lpdf* lpdf; /**< log-likelihood component */
  gsl_vector* e; /**< observed residuals*/

  mcmclib_iwishart_lpdf* w;
} mcmclib_mcar_model;


/** Alloc a new MCAR model object */
mcmclib_mcar_model* mcmclib_mcar_model_alloc(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e);
/** De-alloc a MCAR model object */
void mcmclib_mcar_model_free(mcmclib_mcar_model* p);

/** alpha12sigma full conditional log-distribution. Domain: real line */
double mcmclib_mcar_model_alpha12sigma_lpdf(void* in_p, gsl_vector* alpha12sigma);
/** Gamma full conditional log-distribution.
    Gamma is parametrized as a p x (p-1)/2 + p vector of real valued
    coefficients */
double mcmclib_mcar_model_alphasigma_lpdf(void* in_p, gsl_vector* alphasigma);

/**@}*/
/**@}*/
#endif
