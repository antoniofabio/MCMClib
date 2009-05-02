/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MH_H__
#define __MH_H__

/**\addtogroup metropolis_samplers
@{
\defgroup MH Metropolis-Hastings sampling
@{*/

#include "common.h"
#include "mh_q.h"

/** pointer to a distribution function */
typedef double (*distrfun_p) (void* data, gsl_vector* x);

/**\brief Generic Metropolis-Hastings sampler */
typedef struct {
  gsl_rng* r; /**< rng*/
  distrfun_p logdistr; /**< target log-density fun*/
  void* logdistr_data; /**< target log-density data*/
  mcmclib_mh_q* q; /**< proposal kernel object*/
  gsl_vector* x; /**< current chain value*/
  gsl_vector* x_old; /**< old chain value*/
  double logdistr_old; /**< log-distrib. of old chain value */
  int flag; /**< do we have started sampling? */
  int last_accepted; /**< flag: last move has been accepted?*/
} mcmclib_mh;

/**\brief Alloc a new Metropolis-Hastings sampler
@param r RNG
@param logdistr target prob. density log-function
@param logdistr_data target prob. density fun. extra data
@param q proposal kernel
@param x current chain value
@return the new allocated MH object
*/
mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r,
			     distrfun_p logdistr, void* logdistr_data,
			     mcmclib_mh_q* q, gsl_vector* x);

/**\brief Free a previously allocated MH object*/
void mcmclib_mh_free(mcmclib_mh* p);

/**\brief update chain value
\return 1 if move accepted, 0 if rejected*/
int mcmclib_mh_update(mcmclib_mh* p);

/**\brief Generic (non-symmetric) M-H step \internal*/
int mcmclib_mh_generic_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x,
			    distrfun_p logdistr, void* data,
			    double* plogdistr_old, mcmclib_mh_q* q);
/**@}*/
/**@}*/

#endif
