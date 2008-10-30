#ifndef __METROPOLIS_H__
#define __METROPOLIS_H__
#include "common.h"

/*data: extra data
  x: conditioning value
  y: new value*/
typedef double (*proposal_distr_t) (void* data, gsl_vector* x, gsl_vector* y);

/** Some common M-H support functions */

/**returns 1 for accept, 0 for reject
@param r GSL RNG
@param old old value
@param x vector holding current value (will be eventually updated!)
@param logdistr ptr to log-distribution function
@param data extra data for 'logdistr'
*/
int mcmclib_metropolis_symmetric_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x, distrfun_p logdistr, void* data);

/**returns 1 for accept, 0 for reject
@param r GSL RNG
@param old old value
@param x vector holding current value (will be eventually updated!)
@param logdistr ptr to log-distribution function
@param data extra data for 'logdistr'
@param q proposal log-distribution
@param q_data extra data for 'q'
*/
int mcmclib_metropolis_generic_step(const gsl_rng* r, gsl_vector* old,
				    gsl_vector* x, distrfun_p logdistr, void* data,
				    proposal_distr_t q, void* q_data);

#endif
