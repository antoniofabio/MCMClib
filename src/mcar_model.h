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

/**\brief MCAR(\tilde B, \Gamma) Model

TODO
*/
typedef struct {
  mcmclib_mcar_tilde_lpdf* lpdf;

  gsl_vector* e; /**< observed residuals*/
} mcmclib_mcar_model;

mcmclib_mcar_model* mcmclib_mcar_model_alloc(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e);
void mcmclib_mcar_model_free(mcmclib_mcar_model* p);

double mcmclib_mcar_model_alpha1_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha1);
double mcmclib_mcar_model_alpha2_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha2);
double mcmclib_mcar_model_sigma_lpdf(mcmclib_mcar_model* p, gsl_vector* sigma);
double mcmclib_mcar_model_Gamma_lpdf(mcmclib_mcar_model* p, gsl_vector* gamma);

/**@}*/
/**@}*/
#endif
