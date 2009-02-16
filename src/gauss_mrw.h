#ifndef __GAUSS_MRW_H__
#define __GAUSS_MRW_H__

/**\addtogroup metropolis_samplers
@{
\defgroup GAUSS_MRW Multivariate Gaussian Random Walk
*/
#include "common.h"
#include "mh.h"

/**\brief Multivariate Gaussian Random Walk*/
typedef struct {
  mcmclib_mh* mh;

  gsl_matrix* sigma_prop; /**< proposal covariance matrix*/
} mcmclib_gauss_mrw;

/** alloc (and init) extra Gaussian RW data
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param start_x starting value
@param data extra data to be passed to the distribution function
@param sigma_prop gaussian proposal covariance matrix
*/
mcmclib_gauss_mrw* mcmclib_gauss_mrw_alloc(gsl_rng* r,
	distrfun_p logdistr, void* data, gsl_vector* start_x, const gsl_matrix* sigma_prop);

/** free extra Gaussian RW data*/
void mcmclib_gauss_mrw_free(mcmclib_gauss_mrw* p);

/** Gaussian random walk step*/
int mcmclib_gauss_mrw_update(mcmclib_gauss_mrw* p);

/** GRW proposal log-density (fake) \internal*/
double mcmclib_gauss_mrw_qd(void* ignore, gsl_vector* x, gsl_vector* y);

/** GRW proposal sampler \internal
@param in_p ptr to an mcmclib_gauss_mrw object
@param x object in which to put result
*/
void mcmclib_gauss_mrw_sample(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
