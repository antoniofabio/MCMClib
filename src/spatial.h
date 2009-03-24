#ifndef __MCMCLIB_SPATIAL_H__
#define __MCMCLIB_SPATIAL_H__

#include "mvnorm.h"
/**\addtogroup distributions
 @{*/

/**\addtogroup spatial
 @{*/

/**\brief Multivariate Gaussian distribution with spatial correlation */
typedef struct {
  gsl_vector* mu; /**< mean*/
  gsl_matrix* D; /**< pairwise distances between points*/
  gsl_vector* rho; /**< range*/
  gsl_vector* sigma; /**< sill*/
  gsl_vector* tausq; /**< nugget*/

  gsl_matrix* Sigma;
  mcmclib_mvnorm_lpdf* norm;
} mcmclib_spatial_lpdf;

/** Alloc extra data for a spatial distribution
@param mu mean
@param rho range
@param sigma sill
@param tausq nugget
@param D pairwise distances
*/
mcmclib_spatial_lpdf* mcmclib_spatial_lpdf_alloc(gsl_vector* mu,
						 gsl_vector* rho,
						 gsl_vector* sigma,
						 gsl_vector* tausq,
						 gsl_matrix* D);

/** Free extra data for a spatial distribution
@param p pointer to distrib extra data
*/
void mcmclib_spatial_lpdf_free(mcmclib_spatial_lpdf* p);

/** Spatial log-distribution
@param in_p extra data, allocated via \ref mcmclib_spatial_lpdf_alloc
@return log-pdf
*/
double mcmclib_spatial_lpdf_compute(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
