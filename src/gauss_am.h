#ifndef __GAUSS_AM_H__
#define __GAUSS_AM_H__

/**\addtogroup adaptive
@{
\defgroup GAUSS_AM Gaussian random walk Adaptive Metropolis
@{*/

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "amh.h"
#include "gauss_mrw.h"

/**\brief Adaptive Metropolis Gaussian random walk*/
typedef struct {
  mcmclib_amh* amh;

  /*AM specific fields*/
  mcmclib_gauss_mrw* mrw;
  gsl_matrix* sigma_zero;
  int t0;
  gsl_vector* mean;
  gsl_matrix* cov;
  double sf;
} mcmclib_gauss_am;

/** alloc (and init) extra AM data
@param r RNG state
@param logdistr pointer to a log-distribution function
@param logdistrib_data extra data to be passed to the log-distribution function
@param start_x starting value
@param sigma_zero starting proposal covariance matrix
@param t0 burn-in length before starting adaptation
*/
mcmclib_gauss_am* mcmclib_gauss_am_alloc(gsl_rng* r,
					 distrfun_p logdistr, void* logdistr_data,
					 gsl_vector* start_x,
					 const gsl_matrix* sigma_zero, int t0);
/** free extra AM data*/
void mcmclib_gauss_am_free(mcmclib_gauss_am* p);

/** Adaptive Metropolis Gaussian random walk*/
int mcmclib_gauss_am_update(mcmclib_gauss_am* p);

/** AM gamma update function \internal*/
void mcmclib_gauss_am_update_gamma(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
