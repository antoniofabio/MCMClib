#ifndef __MH_H__
#define __MH_H__

/**\addtogroup metropolis_samplers
@{
\defgroup MH Metropolis-Hastings sampling
*/

#include "common.h"
#include "mh_q.h"

/**\brief Generic Metropolis-Hastings sampler */
typedef struct {
  gsl_rng* r; /**< rng*/
  distrfun_p logdistr; /**< target log-density fun*/
  void* logdistr_data; /**< target log-density data*/
  mcmclib_mh_q* q; /**< proposal kernel object*/
  gsl_vector* x; /**< current chain value*/
  gsl_vector* x_old; /**< old chain value*/
  int last_accepted; /**< flag: last move has been accepted?*/
} mcmclib_mh;

mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r,
			     distrfun_p logdistr, void* logdistr_data,
			     mcmclib_mh_q* q, gsl_vector* x);

void mcmclib_mh_free(mcmclib_mh* p);

/**\brief update chain value
\return 1 if move accepted, 0 if rejected*/
int mcmclib_mh_update(mcmclib_mh* p);

/**@}*/
/**@}*/

#endif
