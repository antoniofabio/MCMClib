#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include "mh.h"

/**\addtogroup metropolis_samplers
@{
\defgroup gauss_rw gauss_rw
@{*/

typedef struct {
  gsl_rng* r;
  double step_size;
} mcmclib_gauss_rw_gamma;

mcmclib_gauss_rw_gamma* mcmclib_gauss_rw_gamma_alloc(gsl_rng* r, double step_size);
void mcmclib_gauss_rw_gamma_free(mcmclib_gauss_rw_gamma* p);
void mcmclib_gauss_rw_sample(void* in_p, gsl_vector* x);
double mcmclib_gauss_rw_qd(void* ignore, gsl_vector* x, gsl_vector* y);

mcmclib_mh_q* mcmclib_gauss_rw_q_alloc(gsl_rng* r, double step_size);
void mcmclib_gauss_rw_q_free(mcmclib_mh_q* p);

/** alloc (and init) Gaussian RW object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param start_x starting value
@param data extra data to be passed to the distribution function
@param step_size gaussian proposal width (s.d.)
*/
mcmclib_mh* mcmclib_gauss_rw_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* data,
				   gsl_vector* start_x, double step_size);

/** free Gaussian RW object*/
void mcmclib_gauss_rw_free(mcmclib_mh* p);

/**@}*/
/**@}*/
#endif
