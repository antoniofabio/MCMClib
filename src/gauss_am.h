#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "common.h"
#include "mvnorm.h"

/** Adaptive Metropolis Gaussian random walk extra data
*/
typedef struct {
	gsl_matrix* sigma_zero;
	int t0;
	gsl_vector* mean;
	gsl_matrix* cov;
	int t;
	gsl_vector* old;
} mcmclib_gauss_am_data;

/** alloc (and init) extra AM data
@param sigma_zero starting proposal covariance matrix
@param t0 burn-in length before starting adaptation
*/
mcmclib_gauss_am_data* mcmclib_gauss_am_alloc(const gsl_matrix* sigma_zero, int t0);
/** free extra AM data
*/
void mcmclib_gauss_am_free(mcmclib_gauss_am_data* p);

/** Adaptive Metropolis Gaussian random walk
@param r RNG state
@param loglik pointer to a log-likelihood function
@param x current point value
@param data extra data to be passed to the log-likelihood function
@param extra pointer to extra AM data
*/
int mcmclib_gauss_am(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	mcmclib_gauss_am_data*);

#endif
