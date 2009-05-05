/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MATRIX_H__
#define __MCMCLIB_MATRIX_H__
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/**\addtogroup misc
 @{*/

/** Compute matrix inverse in-place, by LU factorization */
void mcmclib_matrix_inverse(gsl_matrix* A);

/** check if vector 'v' is sorted in descending order */
int mcmclib_vector_is_sorted_desc(gsl_vector* v);

/** Print vector 'v' on stdout */
void mcmclib_vector_printf(gsl_vector* v);

/** Print matrix 'A' on stdout */
void mcmclib_matrix_printf(gsl_matrix* A);

/**@}*/
#endif
