/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MH_Q_H__
#define __MH_Q_H__

/**\addtogroup metropolis_samplers
@{
\defgroup MH Metropolis-Hastings sampling
@{*/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

typedef void (*free_fun_t) (void* data);

struct mcmclib_mh_q_t;

/** pointer to a proposal log-distribution function
@param data extra data
@param  x conditioning value
@param  y new value
*/
typedef double (*proposal_distr_fun_t) (struct mcmclib_mh_q_t* data,
					gsl_vector* x, gsl_vector* y);

/** pointer to a sampler function
@param data ptr to an mcmclib_mh_q object
*/
typedef void (*sampler_fun_t) (struct mcmclib_mh_q_t* data, gsl_vector* x);

/**\brief Metropolis-Hastings proposal kernel*/
typedef struct mcmclib_mh_q_t {
  gsl_rng* r; /**< RNG used by the sampler function*/
  sampler_fun_t rq; /**< proposal sampler fun*/
  proposal_distr_fun_t dq; /**< proposal density fun*/
  void* gamma; /**< extra proposal kernel data*/
  free_fun_t free_gamma_fun; /**< optional gamma de-allocator fun*/
} mcmclib_mh_q;

/**\brief alloc a new M-H proposal kernel
   @param r RNG
   @param rq proposal sampler fun
   @param dq proposal density fun
   @param gamma misc kernel parameters extra data
   @param free_gamma_fun gamma de-allocator function
 */
mcmclib_mh_q* mcmclib_mh_q_alloc(gsl_rng* r,
				 sampler_fun_t rq,
				 proposal_distr_fun_t dq,
				 void* gamma, free_fun_t free_gamma_fun);
void mcmclib_mh_q_free(mcmclib_mh_q* p);

/**\brief sample a new proposal point*/
void mcmclib_mh_q_sample(mcmclib_mh_q* p, gsl_vector* x);
/**\brief wrapper around the proposal kernel log-density function*/
double mcmclib_mh_q_logd(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y);
/**\brief compute the M-H ratio offset using the kernel q function*/
double mcmclib_mh_q_ratio_offset(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y);

/**@}*/
/**@}*/

#endif
