#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "common.h"

/** Adaptive Metropolis Gaussian random walk
@param r RNG state
@param loglik pointer to a log-likelihood function
@param x current point value
@param data extra data to be passed to the log-likelihood function
@param sigma_zero starting proposal covariance matrix
@param t0 burn-in length before starting adaptation
@param cov pointer to sample covariance matrix observed so far
@param mean pointer to sample mean observed so far
@param t pointer to the number of iterations performed so far
*/
int mcmclib_gauss_am(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	gsl_matrix* sigma_zero, int t0,
	gsl_matrix* cov, int* t);

#endif
