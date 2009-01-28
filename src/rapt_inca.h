#ifndef __RAPT_INCA_H__
#define __RAPT_INCA_H__

#include "rapt.h"

/**\addtogroup adaptive
@{*/
/**\defgroup INCA_RAPT inca_rapt
\brief INCA Regional AdaPTive
@{*/

/** \brief INCA_RAPT sampler data */
typedef struct {
  mcmclib_rapt** sm; /**< array of parallel RAPT samplers*/
  int M; /**< number of parallel chains*/
  gsl_vector** means; /**< array of regions means*/
  gsl_matrix** variances; /**< array of regions variances*/
  gsl_vector* global_mean;
  gsl_matrix* global_variance;
  int t; /**< total number of iterations so far*/
  gsl_vector* n; /**< number of visits in each region*/
} mcmclib_inca_rapt;

/** alloc a new INCA RAPT sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x array of current chain values
@param M number of parallel chains
@param t0 burn-in length before starting adaptation
@param sigma_whole global proposal covariance matrix
@param K number of regions
@param sigma_local array of local proposal covariance matrices
@param which_region boundary computing function
@param which_region_data ptr to extra 'which_region' data
*/
mcmclib_inca_rapt* mcmclib_inca_rapt_alloc(gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 gsl_vector** x, int M, int t0,
				 const gsl_matrix* sigma_whole,
				 int K, gsl_matrix** sigma_local,
				 region_fun_t which_region,
				 void* which_region_data);

/** free  data*/
void mcmclib_inca_rapt_free(mcmclib_inca_rapt* p);

/** Update current value of the INCA_RAPT chains*/
int mcmclib_inca_rapt_update(mcmclib_inca_rapt* p);

/**@}*/
/**@}*/
#endif
