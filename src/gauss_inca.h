#ifndef __GAUSS_INCA_H__
#define __GAUSS_INCA_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "common.h"
#include "mvnorm.h"

/** INCA chains shared data structure
*/
typedef struct {
	gsl_matrix* Sigma_zero;
	int t0;
	int K; /**number of parallel chains*/
	gsl_vector** mean;
	gsl_matrix** variance;
	int* t; /**number of iterates for each chain (on which mean and variance h.b. computed*/
	gsl_vector* mean_global;
	gsl_matrix* variance_global;
	int id; /**current chain id*/
	double sf;
	gsl_matrix* sigma_proposal;
} mcmclib_gauss_inca_pool;

/** INCA chains shared data structure allocator
@param Sigma_zero starting covariance guess
@param t0 burn-in length
@param K number of parallel chains to account for
*/
mcmclib_gauss_inca_pool* mcmclib_gauss_inca_pool_alloc(gsl_matrix* Sigma_zero,int t0,int K);
/** update global mean and variance infos of INCA pool data structure */
void mcmclib_gauss_inca_pool_update_variance(mcmclib_gauss_inca_pool* p);
/** free INCA chains shared memory space */
void mcmclib_gauss_inca_pool_free(mcmclib_gauss_inca_pool* p);

/** INter-Chain Adaptive Gaussian random walk extra data
*/
typedef struct {
	mcmclib_gauss_inca_pool* p;
	gsl_vector* old;
	int id;
} mcmclib_gauss_inca;

/** alloc (and init) extra INCA data
@param pool already allocated INCA shared data
*/
mcmclib_gauss_inca* mcmclib_gauss_inca_alloc(mcmclib_gauss_inca_pool* pool);
/** free extra AM data
*/
void mcmclib_gauss_inca_free(mcmclib_gauss_inca* p);

/** Adaptive Metropolis Gaussian random walk
@param extra pointer to extra AM data
@param r RNG state
@param logdistr pointer to a log-distribution function
@param x current point value
@param data extra data to be passed to the log-distribution function
@note You're supposed to call all the chains referring to the same pool in sequence.
Update of global mean and variance/covariance values automatically occurs upon calling
this function on the last chain!
*/
int mcmclib_gauss_inca_update(mcmclib_gauss_inca* extra, const gsl_rng* r,
	distrfun_p logdistr, gsl_vector* x, void* data);

#endif
