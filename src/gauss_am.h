#ifndef __GAUSS_AM_H__
#define __GAUSS_AM_H__

/**\addtogroup adaptive
@{
\defgroup GAUSS_AM Gaussian random walk Adaptive Metropolis
@{*/

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "amh.h"
#include "gauss_mrw.h"

/**\brief Gaussian AM cumulated sufficient statistics and support data*/
typedef struct {
  gsl_vector* sum_x; /**< cumulated sum of xs*/
  gsl_matrix* sum_xx; /**< cumulated sum of xxs*/
  gsl_matrix* Sigma_eps; /**< pos. definiteness cov. correction additive constant*/
  gsl_matrix* Sigma_zero; /**< starting proposal covariance matrix*/
  int t0; /**< burn in before starting adaptation*/
  double sf; /**< scaling factor*/
} mcmclib_gauss_am_suff;

/** alloc (and init) extra AM data
@param r RNG state
@param logdistr pointer to a log-distribution function
@param logdistrib_data extra data to be passed to the log-distribution function
@param start_x starting value
@param sigma_zero starting proposal covariance matrix
@param t0 burn-in length before starting adaptation
*/
mcmclib_amh* mcmclib_gauss_am_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* logdistr_data,
				    gsl_vector* start_x,
				    const gsl_matrix* sigma_zero, int t0);
/** free extra AM data*/
void mcmclib_gauss_am_free(mcmclib_amh* p);

/** AM gamma update function \internal
@param in_p ptr to an mcmclib_amh object
@param x new sampled data point
*/
void mcmclib_gauss_am_update_gamma(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
