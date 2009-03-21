#ifndef __INCA_RAPTOR_H__
#define __INCA_RAPTOR_H__

#include "inca.h"
#include "raptor.h"

/**\addtogroup adaptive
@{*/
/**\defgroup INCA_RAPT inca_raptor
\brief INCA Regional AdaPTive OnLine Recursive
@{*/

/** alloc a new INCA RAPTOR sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x array of current chain values
@param M number of parallel chains
@param t0 burn-in length before starting adaptation
@param Sigma_zero starting variance-covariance matrix
@param beta_hat starting mixture weights estimates
@param mu_hat starting mixture means estimates
@param Sigma_hat starting mixture variance estimates
*/
mcmclib_inca* mcmclib_inca_raptor_alloc(gsl_rng* r,
					distrfun_p logdistr, void* logdistr_data,
					gsl_vector** x, int M,
					int t0, gsl_matrix* Sigma_zero,
					gsl_vector* beta_hat,
					gsl_vector** mu_hat,
					gsl_matrix** Sigma_hat);

/** free  data*/
void mcmclib_inca_raptor_free(mcmclib_inca* p);

/**@}*/
/**@}*/
#endif
