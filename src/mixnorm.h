#ifndef __MCMCLIB_MIXNORM__
#define __MCMCLIB_MIXNORM__

#include <gsl/gsl_vector.h>
#include "mvnorm.h"

typedef struct {
  gsl_vector* w;
  mcmclib_mvnorm_lpdf** pis;
} mcmclib_mixnorm_lpdf;

/** alloc multivariate gaussian mixture distribution
@param w vector of weights
@param pis array of previously allocated log-pdf functions
*/
mcmclib_mixnorm_lpdf* mcmclib_mixnorm_lpdf_alloc(gsl_vector* w, mcmclib_mvnorm_lpdf** pis);

/** free multivariate gaussian mixture distribution
@param p pointer to distrib data
*/
void mcmclib_mixnorm_lpdf_free(mcmclib_mixnorm_lpdf* p);

/** multivariate gassian mixture log-distribution
@param p data allocated via 'mcmclib_mixnorm_lpdf_alloc'
@return log-pdf
*/
double mcmclib_mixnorm_lpdf_compute(void* p, gsl_vector* x);

#endif
