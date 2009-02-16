#ifndef __INCA_H__
#define __INCA_H__

/**\addtogroup adaptive
@{
\defgroup inca INter-Chain Adaptation
*/

#include "amh.h"

/**\brief Generic INCA sampler */
typedef struct {
  mcmclib_amh* amh;
  gsl_vector** x; /**< array of current chain values*/
  int M; /**< number of parallel chains*/
} mcmclib_inca;

/**\brief alloc a new INCA sampler object
@param update_gamma function which updates current gamma value
@param x array of M current chain values
@param M number of parallel chains
*/
mcmclib_inca* mcmclib_inca_alloc(gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 proposal_distr_t qd_fun, void* qd_data,
				 samplerfun_p q_sampler, void* q_sampler_data,
				 mcmclib_amh_update_gamma_p update_gamma,
				 void* update_gamma_data,
				 gsl_vector** x, int M);

void mcmclib_inca_free(mcmclib_inca* p);

/**\brief update all chains values sequentially
\return 1 if last chain move accepted, 0 if rejected*/
int mcmclib_inca_update(mcmclib_inca* p);

/**@}*/
/**@}*/

#endif
