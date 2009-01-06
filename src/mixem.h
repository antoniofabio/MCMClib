#ifndef __MCMCLIB_MIXEM_H__

/**\addtogroup misc
@{*/

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

/**Fitting a gaussian mixture distribution by the EM algorithm
@param X matrix of observed values
@param K number of mixture components
@param mu array of current mixture components means
@param Sigma array of current mixture components variances
@param P probability of belonging to each mixt. component, one row per point
@param w mixture weights
@param NITER desired number of iterations
*/
void mcmclib_mixem_fit(gsl_matrix* X,
		       gsl_vector** mu, gsl_matrix** Sigma,
		       gsl_vector* w, int NITER);

/**@}*/
#endif