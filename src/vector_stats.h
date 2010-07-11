/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __VECTOR_STATS_H__
#define __VECTOR_STATS_H__
#include <gsl/gsl_statistics_double.h>
#include "common.h"

/**\addtogroup misc
@{*/

/** column means */
void mcmclib_matrix_colmeans(const gsl_matrix* m, gsl_vector* out);
/** row means */
void mcmclib_matrix_rowmeans(const gsl_matrix* m, gsl_vector* out);

/** get variance/covariance matrix out of the 'vertical' matrix 'm' */
void mcmclib_matrix_covariance(const gsl_matrix* m, gsl_matrix* out);

/** update covariance value 'recursively'
@param cov current covariance matrix
@param mean current mean
@param n current sample size
@param x new value
@return nothing. 'cov', 'mean' and 'n' values will be updated as a side-effect
*/
void mcmclib_covariance_update(gsl_matrix* cov, gsl_vector* mean, size_t* n, const gsl_vector* x);

/**Pooled weighted variance
@param beta first component's weight
@param means array of 2 means
@param variances array of two variances
@param output result
*/
void mcmclib_pooled_variance(const double beta, const gsl_vector** means,
			     const gsl_matrix** variances, gsl_matrix* V);

/**check if vector contains only finite values*/
int mcmclib_vector_finite(gsl_vector* x);

/**@}*/

#endif
