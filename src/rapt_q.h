#ifndef __RAPT_Q_H__
#define __RAPT_Q_H__

#include "amh.h"
#include "mvnorm.h"

/**\addtogroup adaptive
@{*/
/**\addtogroup RAPT
@{*/

/** Pointer to a region-computing function */
typedef int (*region_fun_t) (gsl_vector*, void*);

/**\brief RAPT proposal kernel parameters*/
typedef struct {
  region_fun_t which_region; /**< region computing function*/
  void* which_region_data; /**< extra data for the region computing function*/
  int K; /**< number of regions*/

  gsl_matrix* sigma_whole; /**< global proposal covariance matrix*/
  gsl_matrix** sigma_local; /**< array of local proposal covariance matrices*/
  gsl_matrix* lambda; /**< K+1 weights for local and global proposals, in each region*/

  int which_proposal; /**< which proposal have been used in last step*/
  int which_region_x, which_region_old; /**< region info*/

  gsl_vector* workspace; /**< internal utility vector*/
  gsl_vector* q_mean; /**< extra data for (mixture) proposal densities comp.*/
  mcmclib_mvnorm_lpdf** q_k;/**< extra data for (mixture) proposal densities comp.*/
} mcmclib_rapt_gamma;

/**< alloc a new RAPT proposal kernel*/
mcmclib_mh_q* mcmclib_rapt_q_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* logdistr_data,
				   const gsl_matrix* sigma_whole,
				   int K,
				   gsl_matrix** sigma_local,
				   region_fun_t which_region,
				   void* which_region_data);

/**< free a previously allocated RAPT kernel*/
void mcmclib_rapt_q_free(mcmclib_mh_q* p);

/** customly set global proposal weight (same for all regions)*/
void mcmclib_rapt_gamma_set_alpha(mcmclib_rapt_gamma* p, double alpha);

/**@}*/
/**@}*/
#endif
