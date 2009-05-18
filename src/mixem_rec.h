/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MIXEM_REC_H__
#define __MCMCLIB_MIXEM_REC_H__
/** \addtogroup misc
@{
\addtogroup mixemrec Recursive mixture fitting
@{
*/

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "mvnorm.h"

/**Gaussians mixture fitting by recursive EM*/
typedef struct {
  /*current parameters estimates*/
  gsl_vector** mu;
  gsl_matrix** Sigma;
  gsl_vector* beta;

  /*extra data*/
  int n;
  mcmclib_mvnorm_lpdf** pi_k;
  gsl_matrix** X_sq_sum;
  gsl_vector** X_sum;
  gsl_vector* beta_sum;
  gsl_vector* beta_i;
} mcmclib_mixem_rec;

/**alloc \ref mcmclib_mixem_rec data.
@param mu array of current means estimates
@param Sigma array of current variances estimates
@param beta vector of current weights estimates
*/
mcmclib_mixem_rec* mcmclib_mixem_rec_alloc(gsl_vector** mu,
					   gsl_matrix** Sigma,
					   gsl_vector* beta);

/**free mixem_rec data*/
void mcmclib_mixem_rec_free(mcmclib_mixem_rec* p);

/**accumulate a new datapoint infos to mixem_rec data*/
void mcmclib_mixem_rec_add(mcmclib_mixem_rec* p, gsl_vector* y);

/**update mixem_rec estimates using currently accumulated data
   @return GSL_SUCCESS if all goes right
 */
int mcmclib_mixem_rec_update(mcmclib_mixem_rec* p);

/**@}@}*/
#endif
