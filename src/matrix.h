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

/** add and scale matrices

    dest = alpha * (A + B) */
void mcmclib_matrix_addscale(gsl_matrix* dest,
			     const gsl_matrix* A, const gsl_matrix* B, double alpha);

/** Compute Cholesky dec. in place.
    @return !GSL_SUCCESS if the cholesky decomposition fails */
int mcmclib_cholesky_decomp(gsl_matrix* A);

/** Compute matrix inverse in-place, by Cholesky dec.
    @return !GSL_SUCCESS if the cholesky decomposition fails */
int mcmclib_cholesky_inverse(gsl_matrix* A);

/** Trace of the logarithm of matrix 'A' */
double mcmclib_matrix_logtrace(const gsl_matrix* A);

/** Compute eigenvalues of a real, symmetric matrix*/
void mcmclib_matrix_symm_eigenvalues(const gsl_matrix* A, gsl_vector* out);

/** Compute matrix inverse in-place, by LU factorization */
void mcmclib_matrix_inverse(gsl_matrix* A);

/** check if all elements of 'v' are finite real numbers */
int mcmclib_vector_is_finite(const gsl_vector* x);

/** check if vector 'v' is sorted in descending order */
int mcmclib_vector_is_sorted_desc(gsl_vector* v);

/** Print vector 'v' on stdout */
void mcmclib_vector_printf(gsl_vector* v);

/** Print matrix 'A' on stdout */
void mcmclib_matrix_printf(gsl_matrix* A);

/**@}*/
#endif
