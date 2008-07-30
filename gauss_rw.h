#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#define isfinite(x) !(isnan((x)) || isinf((x)))

/** Gaussian random walk
@param r RNG state
@param loglik pointer to a log-likelihood function
@param x current point value
@param data extra data to be passed to the log-likelihood function
*/
int mcmclib_gauss_rw(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	const double step_size);

#endif
