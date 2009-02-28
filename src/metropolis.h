#ifndef __METROPOLIS_H__
#define __METROPOLIS_H__
#include "common.h"
#include "mh_q.h"

/**\addtogroup metropolis_samplers
@{
*/

/** M-H step for a symmetric proposal density
@param r GSL RNG
@param old old value
@param x vector holding current value (will be eventually updated!)
@param logdistr ptr to log-distribution function
@param data extra data for 'logdistr'
@return 1 if accept, 0 if reject
*/
int mcmclib_metropolis_symmetric_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x, distrfun_p logdistr, void* data);

/** M-H step for a generic (non-symmetric) proposal density
@param r GSL RNG
@param old old value
@param x vector holding current value (will be eventually updated!)
@param logdistr ptr to log-distribution function
@param data extra data for 'logdistr'
@param q proposal log-distribution
@param q_data extra data for 'q'
@return 1 if accept, 0 if reject
*/
int mcmclib_metropolis_generic_step(const gsl_rng* r, gsl_vector* old,
				    gsl_vector* x, distrfun_p logdistr, void* data,
				    proposal_distr_t q, void* q_data);

/**@}*/
#endif
