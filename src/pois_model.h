/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_POIS_MODEL_H__
#define __MCMCLIB_POIS_MODEL_H__
/**\addtogroup Models
 @{*/
/**\addtogroup Poisson
 @{*/

#include <gsl/gsl_matrix.h>
#include "amh.h"

typedef struct {
  gsl_vector* beta;

  /*internal*/
  gsl_matrix* X; /**< regression matrix */
  const gsl_vector* y; /**< observed values */
  const gsl_vector* offset; /**< (optional) offset on the mean */
  gsl_vector* b0; /**< prior beta mean */
  gsl_matrix* B0; /**< prior beta precision */
  gsl_vector* mu; /**< log-means (= X beta) */
  double ldet; /**< log-determinant of the inverse of B0 */
  gsl_vector *work1, *work2; /**< workspace memory */
} mcmclib_pois_model;

mcmclib_pois_model* mcmclib_pois_model_alloc(const gsl_matrix* X, const gsl_vector* y);
void mcmclib_pois_model_free(mcmclib_pois_model* p);
int mcmclib_pois_model_set_prior_mean(mcmclib_pois_model* p, const gsl_vector* b0);
int mcmclib_pois_model_set_prior_var(mcmclib_pois_model* p, const gsl_matrix* B0);
int mcmclib_pois_model_set_offset(mcmclib_pois_model* p, const gsl_vector* offset);

double mcmclib_pois_model_llik(mcmclib_pois_model* p, gsl_vector* x);
double mcmclib_pois_model_lprior(mcmclib_pois_model* p, gsl_vector* x);
double mcmclib_pois_model_lpdf(void* in_p, gsl_vector* x);

typedef struct {
  mcmclib_pois_model* model;
  mcmclib_amh* sampler;
} mcmclib_pmodel_sampler;

mcmclib_pmodel_sampler* mcmclib_pmodel_sampler_alloc(const gsl_matrix* X,
						     const gsl_vector* y,
						     const gsl_vector* offset,
						     gsl_rng* rng,
						     double sigma0,
						     int burnin);

void mcmclib_pmodel_sampler_free(mcmclib_pmodel_sampler* p);

int mcmclib_pmodel_sampler_update(mcmclib_pmodel_sampler* p);

/**@}*/
/**@}*/
#endif
