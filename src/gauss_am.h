#ifndef __GAUSS_AM_H__
#define __GAUSS_AM_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "common.h"
#include "mvnorm.h"

/** Adaptive Metropolis Gaussian random walk extra data
*/
typedef struct {
	/**common MCMC fields*/
	gsl_rng* r;
	distrfun_p logdistr;
	void* logdistr_data;
	gsl_vector* current_x;
	gsl_vector* old;
	/**AM specific fields*/
	gsl_matrix* sigma_zero;
	int t0;
	gsl_vector* mean;
	gsl_matrix* cov;
	int t;
	double sf;
	gsl_matrix* sigma_prop;
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
	distrfun_p logdistr, void* logdistr_data, gsl_vector* start_x,
	const gsl_matrix* sigma_zero, int t0);
/** free extra AM data
*/
void mcmclib_gauss_am_free(mcmclib_gauss_am* p);

/** Adaptive Metropolis Gaussian random walk
@param p pointer to an AM object
*/
int mcmclib_gauss_am_update(mcmclib_gauss_am* p);

#endif
