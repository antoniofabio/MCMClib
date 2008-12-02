#ifndef __MCMC_HIERARCHICAL_H__
#define __MCMC_HIERARCHICAL_H__

#include "common.h"

/** \addtogroup distributions
@{*/

/**\brief extra data for log-posterior distribution function*/
typedef struct {
  gsl_vector* x;
  distrfun_p prior;
  void* parms;
  distrfun_p loglik;
  gsl_vector** childs;
  void** child_parms;
  gsl_vector* workspace;
} mcmclib_post_lpdf;

/** alloc logposterior distrib. extra data
@param x pointer to vector node 'current' value
@param prior	prior distribution
@param parms	prior distr parameters
@param loglik (common) likelihood function
@param childs null-terminated array of pointers to child vectors
@param child_parms vector of parameters for each child vector
*/
mcmclib_post_lpdf* mcmclib_post_lpdf_alloc(gsl_vector* x, distrfun_p prior, void* parms,
	distrfun_p loglik, gsl_vector** childs, void** child_parms);

/** free logposterior extra data*/
void mcmclib_post_lpdf_free(mcmclib_post_lpdf* p);

/** log-posterior distribution function
@param data pointer to a \ref mcmclib_post_lpdf object allocated with
   \ref mcmclib_post_lpdf_alloc*/
double mcmclib_post_lpdf_compute(void* data, gsl_vector* x);

/**@}*/
#endif
