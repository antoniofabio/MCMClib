#ifndef __AMH_H__
#define __AMH_H__

/**\addtogroup adaptive
@{
\defgroup amh Adaptive Metropolis-Hastings sampling
@{*/

#include "mh.h"

/**\brief Proposal kernel param updating function
@param p ptr to a 'suff' data object
@param x newly sampled point
*/
typedef void (*mcmclib_amh_update_gamma_p) (void* p, gsl_vector* x);

/**\brief Generic Adaptive Metropolis-Hastings sampler */
typedef struct {
  mcmclib_mh* mh;
  void* suff; /**< sufficient data accumulated up to current iteration*/
  mcmclib_amh_update_gamma_p update_gamma;
  int n; /**< current iteration number*/
} mcmclib_amh;

/**\brief alloc a new generic AMH sampler
@param mh base MH sampler to be extended
@param suff extra data to be feeded to the update-gamma function
@param update_gamma update proposal kernel parameters
*/
mcmclib_amh* mcmclib_amh_alloc(mcmclib_mh* mh, void* suff,
			       mcmclib_amh_update_gamma_p update_gamma);

/**\brief free previously allocated AMH sampler*/
void mcmclib_amh_free(mcmclib_amh* p);

/**\brief update chain value
\return 1 if move accepted, 0 if rejected*/
int mcmclib_amh_update(mcmclib_amh* p);

/**@}*/
/**@}*/

#endif
