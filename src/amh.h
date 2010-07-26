/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __AMH_H__
#define __AMH_H__

/**\addtogroup adaptive
@{
\defgroup amh Adaptive Metropolis-Hastings sampling
@{*/

#include "mh.h"

struct mcmclib_amh_t;

/**\brief Proposal kernel param updating function */
typedef void (*mcmclib_amh_update_gamma_t) (void* suff, void* gamma);

/**\brief Proposal kernel param updating function */
typedef void (*mcmclib_amh_update_suff_t) (void* suff, const gsl_vector* x);

/**\brief Generic Adaptive Metropolis-Hastings sampler */
typedef struct mcmclib_amh_t {
  mcmclib_mh* mh;
  mcmclib_amh_update_suff_t update_suff; /**< update sufficient statistics*/
  mcmclib_amh_update_gamma_t update_gamma; /**< update proposal kernel params*/
  size_t n; /**< current iteration number*/
  size_t n0; /**< adaptation burnin*/
  void* suff; /**< sufficient data accumulated up to current iteration*/
  free_fun_t free_suff; /**< sufficient data destructor*/
} mcmclib_amh;

/**\brief alloc a new generic AMH sampler
@param mh base MH sampler to be extended
@param n0 num. iterations to accumulate before starting to update the proposal kernel
@param suff suff. stats. data accumulated by the update_suff function
@param free_suff sufficient data destructor
@param update_suff update sufficient stats
@param update_gamma update proposal kernel parameters
*/
mcmclib_amh* mcmclib_amh_alloc(mcmclib_mh* mh, const size_t n0,
			       void* suff,
			       free_fun_t free_suff,
			       mcmclib_amh_update_suff_t update_suff,
			       mcmclib_amh_update_gamma_t update_gamma);

/**\brief free previously allocated AMH sampler*/
void mcmclib_amh_free(mcmclib_amh* p);

/**\brief update chain value
\return 1 if move accepted, 0 if rejected*/
int mcmclib_amh_update(mcmclib_amh* p);

/**\brief reset chain */
void mcmclib_amh_reset(mcmclib_amh* p);

/**@}*/
/**@}*/

#endif
