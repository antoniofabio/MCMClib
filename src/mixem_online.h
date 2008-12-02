#ifndef __MCMCLIB_MIXEM_ONLINE_H__
#define __MCMCLIB_MIXEM_ONLINE_H__

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "mixolem_suff.h"
#include "mvnorm.h"

/** \addtogroup misc Miscellanea
@{
\addtogroup mixolem
@{
*/

/**On-Line EM for gaussian mixtures*/
typedef struct {
  /*pointers to current parameters estimates*/
  mcmclib_mixolem_suff* gamma; /**< pointer to current parameters estimates*/

  mcmclib_mixolem_suff* s; /**< pointer to currently accumulated suff. stat.*/
  mcmclib_mixolem_suff* si; /**< pointer to last point suff. stat.*/

  double eta_eps; /**< learning rate*/

  mcmclib_mvnorm_lpdf** pi_k; /**< array of mixture components densities*/
  int n, n0; /**< current it. number, burn-in length*/

  gsl_vector* beta;
  gsl_vector** mu;
  gsl_matrix** Sigma;
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

/**update mixem_online estimates using the new datapoint*/
void mcmclib_mixem_online_update(mcmclib_mixem_online* p, gsl_vector* y);

/**update sufficient statistic using the new datapoint
\internal*/
void mcmclib_mixem_online_update_s(mcmclib_mixem_online* p, gsl_vector* y);

/**compute gamma estimate basing on a sufficient stat.
\internal*/
void mcmclib_mixem_online_update_gamma(mcmclib_mixolem_suff* gamma,
				       mcmclib_mixolem_suff* s);

/**@}*/
/**@}*/
#endif
