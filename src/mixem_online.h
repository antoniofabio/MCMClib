#ifndef __MCMCLIB_MIXEM_ONLINE_H__
#define __MCMCLIB_MIXEM_ONLINE_H__
/** \file
\brief Fitting a gaussian mixture by an Online EM algorithm
*/

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "mvnorm.h"

typedef struct {
  /*pointers to current parameters estimates*/
  gsl_vector** mu;
  gsl_matrix** Sigma;
  gsl_vector* beta;

  /*pointers to current sufficient statistics estimates*/
  gsl_vector* delta;
  gsl_vector** delta_x;
  gsl_matrix** delta_xx;

  /*learning rate*/
  double eta_eps;

  /*extra data*/
  mcmclib_mvnorm_lpdf** pi_k;
  int n, n0;
  gsl_vector* deltai;
  gsl_vector** delta_xi;
  gsl_matrix** delta_xxi;
} mcmclib_mixem_online;

/**alloc mixem_online data
@param mu array of current means estimates
@param Sigma array of current variances estimates
@param beta vector of current weights estimates
@param eta_eps learning factor is: n^(-0.5 -eta_eps)
@param n0 don't update gamma estimate until n>n0
*/
mcmclib_mixem_online* mcmclib_mixem_online_alloc(gsl_vector** mu,
						 gsl_matrix** Sigma,
						 gsl_vector* beta,
						 double eta_eps,
						 int n0);

/**free mixem_online data*/
void mcmclib_mixem_online_free(mcmclib_mixem_online* p);

/**update mixem_rec estimates using the new datapoint*/
void mcmclib_mixem_online_update(mcmclib_mixem_online* p, gsl_vector* y);

#endif
