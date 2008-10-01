#ifndef __GAUSS_MRW_H__
#define __GAUSS_MRW_H__

#include "common.h"

/** Multivariate Gaussian Random Walk extra data
*/
typedef struct {
	/**common MCMC fields*/
	gsl_rng* r;
	distrfun_p logdistr;
	void* logdistr_data;
	gsl_vector* current_x;
	gsl_vector* old;

	/**MRW specific fields*/
	gsl_matrix* sigma_prop;
} mcmclib_gauss_mrw;

/** alloc (and init) extra Gaussian RW data
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param start_x starting value
@param data extra data to be passed to the distribution function
@param sigma_prop gaussian proposal covariance matrix
*/
mcmclib_gauss_mrw* mcmclib_gauss_rw_alloc(gsl_rng* r,
	distrfun_p logdistr, void* data, gsl_vector* start_x, gsl_matrix* sigma_prop);
/** free extra Gaussian RW data
*/
void mcmclib_gauss_mrw_free(mcmclib_gauss_mrw* p);

/** Gaussian random walk
@param p a MRW object
*/
int mcmclib_gauss_mrw_update(mcmclib_gauss_mrw* p);

#endif
