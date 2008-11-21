#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

/**\file
\brief Gaussian random walk
*/

#include "common.h"

/** Gaussian random walk extra data
*/
typedef struct {
	gsl_rng* r;
	distrfun_p logdistr;
	void* logdistr_data;
	gsl_vector* current_x;
	double step_size;
	/**internal*/
	gsl_vector* old;
} mcmclib_gauss_rw;

/** alloc (and init) Gaussian RW object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param start_x starting value
@param data extra data to be passed to the distribution function
@param step_size gaussian proposal width (s.d.)
*/
mcmclib_gauss_rw* mcmclib_gauss_rw_alloc(gsl_rng* r,
	distrfun_p logdistr, void* data, gsl_vector* start_x, double step_size);

/** free extra Gaussian RW object
*/
void mcmclib_gauss_rw_free(mcmclib_gauss_rw* p);

/** Gaussian random walk step
*/
int mcmclib_gauss_rw_update(mcmclib_gauss_rw* p);

#endif
