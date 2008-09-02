#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include "common.h"

/** Gaussian random walk extra data
*/
typedef struct {
	double step_size;
	gsl_vector* old;
} mcmclib_gauss_rw_data;

/** alloc (and init) extra Gaussian RW data
@param step_size gaussian proposal width (s.d.)
@param dim ambient space dimension
*/
mcmclib_gauss_rw_data* mcmclib_gauss_rw_alloc(double step_size, int dim);
/** free extra Gaussian RW data
*/
void mcmclib_gauss_rw_free(mcmclib_gauss_rw_data* p);

/** Gaussian random walk
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param x current point value
@param data extra data to be passed to the distribution function
*/
int mcmclib_gauss_rw(mcmclib_gauss_rw_data* e, const gsl_rng* r,
	distrfun_p logdistr, gsl_vector* x, void* data);

#endif
