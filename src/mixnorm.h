/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MIXNORM__
#define __MCMCLIB_MIXNORM__

/**\addtogroup distributions
@{
\addtogroup multivariate
@{
*/

#include <gsl/gsl_vector.h>
#include "mvnorm.h"

/** Multivariate gaussian mixture */
typedef struct {
  gsl_vector* w; /**< mixture weights*/
  mcmclib_mvnorm_lpdf** pis; /**< mixture components*/
} mcmclib_mixnorm_lpdf;

/** alloc multivariate gaussian mixture distribution
@param w vector of weights
@param pis array of previously allocated log-pdf functions
@return the new allocated \ref mcmclib_mixnorm_lpdf object
*/
mcmclib_mixnorm_lpdf* mcmclib_mixnorm_lpdf_alloc(gsl_vector* w, mcmclib_mvnorm_lpdf** pis);

/** free multivariate gaussian mixture distribution
@param p pointer to distrib data
*/
void mcmclib_mixnorm_lpdf_free(mcmclib_mixnorm_lpdf* p);

/** multivariate gassian mixture log-distribution
@param p data allocated via 'mcmclib_mixnorm_lpdf_alloc'
@param x point in which to compute the lpdf
@return log-pdf
*/
double mcmclib_mixnorm_lpdf_compute(void* p, gsl_vector* x);

/**
@}
@}
*/
#endif
