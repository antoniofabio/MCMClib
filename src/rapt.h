#ifndef __RAPT_H__
#define __RAPT_H__

#include "common.h"
#include "mvnorm.h"

/**\addtogroup adaptive
@{*/
/**\defgroup RAPT rapt
\brief Regional AdaPTive
@{*/

/** Pointer to a region-computing function */
typedef int (*region_fun_t) (gsl_vector*, void*);

/** \brief Regional Adaptive sampler data */
typedef struct {
  /*common MCMC fields*/
  gsl_rng* r;
  distrfun_p logdistr;
  void* logdistr_data;
  gsl_vector* current_x;
  gsl_vector* old;
  
  /*rapt specific fields*/
  int t0; /**< burn-in length*/
  gsl_matrix* sigma_whole; /**< global proposal covariance matrix*/
  int K; /**< number of regions*/
  gsl_matrix** sigma_local; /**< array of local proposal covariance matrices*/
  region_fun_t which_region; /**< boundary computing function*/
  void* which_region_data; /**< ptr to extra 'which_region' data*/
  gsl_matrix* lambda; /**< K+1 weights for local and global proposals, in each region*/

  /*extra infos*/
  int t; /**<number of iterations done so far*/
  gsl_vector** means; /**< array of regions means*/
  gsl_matrix** variances; /**< array of regions variances*/
  gsl_vector* global_mean;
  gsl_matrix* global_variance;
  gsl_vector* n; /**< number of visits in each region*/
  int which_proposal; /**< which proposal have been used in last step*/
  int accepted; /**< last step was an acceptance or a rejection?*/
  int which_region_x, which_region_old; /**< region info*/

  /*internal data*/
  gsl_matrix* Sigma_eps; /**< additive perturbation factor for variances updating*/
  gsl_vector* workspace; /**< utility workspace memory*/
  gsl_vector* q_mean; /**< extra data for (mixture) proposal densities comp.*/
  mcmclib_mvnorm_lpdf** q_k;/**< extra data for (mixture) proposal densities comp.*/
} mcmclib_rapt;

/** alloc a new RAPT sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x current chain value
@param t0 burn-in length before starting adaptation
@param sigma_whole global proposal covariance matrix
@param K number of regions
@param sigma_local array of local proposal covariance matrices
@param which_region boundary computing function
@param which_region_data ptr to extra 'which_region' data
*/
mcmclib_rapt* mcmclib_rapt_alloc(gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 gsl_vector* x,
				 int t0,
				 const gsl_matrix* sigma_whole,
				 int K,
				 gsl_matrix** sigma_local,
				 region_fun_t which_region,
				 void* which_region_data);

/** free  data*/
void mcmclib_rapt_free(mcmclib_rapt* p);

/** Update current value of a RAPT chain*/
int mcmclib_rapt_update(mcmclib_rapt* p);

/** update local and global proposals covariance matrices
basing on current (region-specific) sample variances*/
void mcmclib_rapt_update_proposals(mcmclib_rapt* p);

/** update local and global proposals covariance matrices
using custom estimated local and global variances
 */
void mcmclib_rapt_update_proposals_custom(mcmclib_rapt* p,
					  gsl_matrix** variances,
					  gsl_matrix* global_variance);

/**@}*/
/**@}*/
#endif
