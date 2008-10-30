#ifndef __RAPT_H__
#define __RAPT_H__

#include "common.h"

typedef int (*region_fun_t) (gsl_vector*, void*);

/** Regional Adaptive sampler data
*/
typedef struct {
  /**common MCMC fields*/
  gsl_rng* r;
  distrfun_p logdistr;
  void* logdistr_data;
  gsl_vector* current_x;
  gsl_vector* old;
  
  /**rapt specific fields*/
  gsl_matrix* sigma_whole; /*global proposal covariance matrix*/
  int K; /*number of regions*/
  gsl_matrix** sigma_local; /*array of local proposal covariance matrices*/
  region_fun_t which_region; /*boundary computing function*/
  void* which_region_data; /*ptr to extra 'which_region' data*/
} mcmclib_rapt;

/** alloc (and init) extra Gaussian RW data
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param data extra data to be passed to the distribution function
@param x current chain value
@param sigma_whole global proposal covariance matrix
@param K number of regions
@param sigma_local array of local proposal covariance matrices
@param which_region boundary computing function
@param which_region_data ptr to extra 'which_region' data
*/
mcmclib_rapt* mcmclib_rapt_alloc(
				 gsl_rng* r,
				 distrfun_p logdistr, void* data,
				 gsl_vector* x,
				 const gsl_matrix* sigma_whole,
				 int K,
				 gsl_matrix** sigma_local,
				 region_fun_t which_region,
				 void* which_region_data);

/** free  data
*/
void mcmclib_rapt_free(mcmclib_rapt* p);

/** Gaussian random walk
@param p a MRW object
*/
int mcmclib_rapt_update(mcmclib_rapt* p);

#endif
