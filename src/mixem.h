#ifndef __MCMCLIB_MIXEM_H__

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

#define MIXEM_MAXITER 10

/**Fitting a gaussian mixture distribution by the EM algorithm
@param X matrix of observed values
@param K number of mixture components
@param mu array of current mixture components means
@param Sigma array of current mixture components variances
@param P probability of belonging to each mixt. component, one row per point
@param w mixture weights
@param eps convergence criteria
*/
void mcmclib_mixem_fit(gsl_matrix* X, int K,
		       gsl_vector** mu, gsl_matrix** Sigma,
		       gsl_matrix* P, gsl_vector* w, double eps);

#endif
