#ifndef __MH_H__
#define __MH_H__

/**\addtogroup metropolis_samplers
@{
\defgroup MH Metropolis-Hastings sampling
*/

#include "common.h"
#include "metropolis.h"

/**\brief Generic Metropolis-Hastings sampler */
typedef struct {
  gsl_rng* r; /**< rng*/
  distrfun_p logdistr; /**< target log-density fun*/
  void* logdistr_data; /**< target log-density data*/
  gsl_vector* x; /**< current chain value*/
  gsl_vector* x_old; /**< old chain value*/
  proposal_distr_t qd_fun; /**< proposal kernel density fun*/
  void* qd_data; /**< proposal kernel density data*/
  distrfun_p q_sampler; /**< proposal kernel fun*/
  void* q_sampler_data; /**< proposal kernel fun data*/
  int last_accepted; /**< flag: last move has been accepted?*/
} mcmclib_mh;

mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r,
			     distrfun_p logdistr,
			     void* logdistr_data,
			     gsl_vector* x,
			     proposal_distr_t qd_fun,
			     void* qd_data,
			     distrfun_p q_sampler,
			     void* q_sampler_data);

void mcmclib_mh_alloc(mcmclib_mh* p);

/**\brief update chain value
\return 1 if move accepted, 0 if rejected*/
int mcmclib_mh_update(mcmclib_mh* p);

/**@}*/
/**@}*/

#endif
