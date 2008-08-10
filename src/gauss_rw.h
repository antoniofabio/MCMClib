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
@param loglik pointer to a log-likelihood function
@param x current point value
@param data extra data to be passed to the log-likelihood function
@param step_size gaussian proposal width (s.d.)
*/
int mcmclib_gauss_rw(const gsl_rng* r,
	distrfun_p loglik, gsl_vector* x, const void* data,
	mcmclib_gauss_rw_data* e);

#endif
