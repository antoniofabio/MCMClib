#ifndef __MCMC_HIERARCHICAL_H__
#define __MCMC_HIERARCHICAL_H__

#include "common.h"

/** extra data for log-posterior function
*/
typedef struct {
	distrfun_p prior;
	void* parms;
	distrfun_p loglik;
	gsl_vector** childs;
	void** child_parms;
} post_lpdf_p;

/** alloc logposterior extra data
@param prior	prior distribution
@param parms	prior distr parameters
@param loglik (common) likelihood function
@param childs null-terminated array of pointers to child vectors
@param child_parms vector of parameters for each child vector
*/
post_lpdf_p* mcmclib_lpdf_post_alloc(distrfun_p prior, void* parms,
	distrfun_p loglik, gsl_vector** childs, void** child_parms);

/** free logposterior extra data
*/
void mcmclib_lpdf_post_free(post_lpdf_p* p);

/** log-posterior distribution function
*/
double mcmclib_lpdf_post(gsl_vector* x, void* data);

#endif
