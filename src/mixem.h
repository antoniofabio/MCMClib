/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MIXEM_H__

/**\addtogroup misc
@{
*/

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

/**Fitting a gaussian mixture distribution by the EM algorithm
@param X matrix of observed values
@param mu array of current mixture components means
@param Sigma array of current mixture components variances
@param w mixture weights
@param NITER desired number of iterations
@return GSL_SUCCESS if all goes right, a GSL error code otherwise
*/
int mcmclib_mixem_fit(gsl_matrix* X,
		      gsl_vector** mu, gsl_matrix** Sigma,
		      gsl_vector* w, size_t NITER);

/**@}*/
#endif
