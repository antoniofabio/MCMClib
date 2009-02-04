#ifndef __AM_INCA_H__
#define __AM_INCA_H__

#include "gauss_mrw.h"

/**\addtogroup adaptive
@{*/
/**\defgroup AM_INCA am_inca
\brief INCA Metropolis
@{*/

/** \brief AM_INCA sampler data */
typedef struct {
  mcmclib_gauss_mrw** sm; /**< array of parallel MH samplers*/
  int M; /**< number of parallel chains*/
  gsl_vector* global_mean;
  gsl_matrix* global_variance;
  int t0; /**< burn-in length*/
  int t; /**< total number of iterations so far*/
  gsl_matrix* Sigma_eps; /**< variance correction factor*/
} mcmclib_am_inca;

/** alloc a new AM INCA sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param start_x starting value
@param M number of parallel chains
@param t0 burn-in length
@param data extra data to be passed to the distribution function
@param sigma_prop starting gaussian proposal covariance matrix
*/
mcmclib_am_inca* mcmclib_am_inca_alloc(gsl_rng* r,
				       distrfun_p logdistr, void* data,
				       gsl_vector** x, int M, int t0,
				       const gsl_matrix* sigma_prop);

/** free  data*/
void mcmclib_am_inca_free(mcmclib_am_inca* p);

/** Update current value of the INCA_RAPT chains*/
int mcmclib_am_inca_update(mcmclib_am_inca* p);

/**@}*/
/**@}*/
#endif