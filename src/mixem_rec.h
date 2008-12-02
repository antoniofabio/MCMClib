#ifndef __MCMCLIB_MIXEM_REC_H__
#define __MCMCLIB_MIXEM_REC_H__
/** \addtogroup misc
@{
\addtogroup mixemrec Recursive mixture fitting (deprecated)
@{
*/

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "mvnorm.h"

/**Gaussians mixture fitting by recursive EM*/
typedef struct {
  /*current parameters estimates*/
  gsl_vector** mu;
  gsl_matrix** Sigma;
  gsl_vector* beta;

  /*extra data*/
  int n;
  mcmclib_mvnorm_lpdf** pi_k;
  gsl_matrix** X_sq_sum;
  gsl_vector** X_sum;
  gsl_vector* beta_sum;
  gsl_vector* beta_i;
} mcmclib_mixem_rec;

/**alloc \ref mcmclib_mixem_rec data.
@param mu array of current means estimates
@param Sigma array of current variances estimates
@param beta vector of current weights estimates
*/
mcmclib_mixem_rec* mcmclib_mixem_rec_alloc(gsl_vector** mu,
					   gsl_matrix** Sigma,
					   gsl_vector* beta);

/**free mixem_rec data*/
void mcmclib_mixem_rec_free(mcmclib_mixem_rec* p);

/**accumulate a new datapoint infos to mixem_rec data*/
void mcmclib_mixem_rec_add(mcmclib_mixem_rec* p, gsl_vector* y);

/**update mixem_rec estimates using currently accumulated data*/
void mcmclib_mixem_rec_update(mcmclib_mixem_rec* p);

/**@}@}*/
#endif
