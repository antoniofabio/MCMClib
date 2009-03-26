#ifndef __MCMCLIB_AR1_LPDF_H__
#define __MCMCLIB_AR1_LPDF_H__

#include "mvnorm.h"
/**\addtogroup distributions
 @{*/

/**\addtogroup AR1
 @{*/

/**\brief Multivariate Gaussian distribution with AR(1) correlations */
typedef struct {
  gsl_vector* mu;
  gsl_vector* sigma;
  gsl_vector* phi;
  mcmclib_mvnorm_lpdf* norm;
} mcmclib_ar1_lpdf;

/** Alloc extra data for an AR(1) distribution
@param mu mean
@param sigma variance
@param phi autocorrelation
*/
mcmclib_ar1_lpdf* mcmclib_spatial_lpdf_alloc(gsl_vector* mu,
					     gsl_vector* sigma,
					     gsl_vector* phi);

/** Free extra data for an AR(1) distribution
*/
void mcmclib_ar1_lpdf_free(mcmclib_ar1_lpdf* p);

/** AR(1) log-distribution
@param in_p extra data, allocated via \ref mcmclib_ar1_lpdf_alloc
@return log-pdf
*/
double mcmclib_ar1_lpdf_compute(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
