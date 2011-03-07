/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MIXOLEM_SUFF_H__
#define __MCMCLIB_MIXOLEM_SUFF_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/** \addtogroup misc
@{
\addtogroup mixolem Gaussian mixture fitting by On-line EM
@{
*/


/** Mixture of multivariate normals complete data sufficient statistics.

In the fields description, \a K is the number of components,
\a d is the ambient space dimension.
*/
typedef struct {
  gsl_vector* delta; /**< vector of K scalars*/
  gsl_vector** delta_x; /**< array of K d-dimensional vectors*/
  gsl_matrix** delta_xx; /**< array of K dxd matrices*/
} mcmclib_mixolem_suff;

/**alloc a new mixolem_suff struct, pointing to externally owned data.

\note You have to free resources allocated by this function
by using the standard \a free*/
mcmclib_mixolem_suff* mcmclib_mixolem_suff_makeref(gsl_vector* delta,
						   gsl_vector** delta_x,
						   gsl_matrix** delta_xx);
/**alloc a new mixolem_suff object*/
mcmclib_mixolem_suff* mcmclib_mixolem_suff_alloc(size_t K, size_t dim);
/**free mixolem_suff object allocated with \ref mcmclib_mixolem_suff_alloc*/
void mcmclib_mixolem_suff_free(mcmclib_mixolem_suff* p);
/**copy contents of \a src in \a dest*/
void mcmclib_mixolem_suff_memcpy(mcmclib_mixolem_suff* dest, mcmclib_mixolem_suff* src);
/**adds \a src to \a dest , put result in \a dest*/
void mcmclib_mixolem_suff_add(mcmclib_mixolem_suff* dest, mcmclib_mixolem_suff* src);
/**scale all \a p contents by \a alpha*/
void mcmclib_mixolem_suff_scale(mcmclib_mixolem_suff* p, double alpha);

/**
@}
@}
*/
#endif
