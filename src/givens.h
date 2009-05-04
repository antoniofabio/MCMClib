/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_GIVENS_H__
#define __MCMCLIB_GIVENS_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/**\addtogroup misc
 @{*/

/** Build an ortogonal matrix from its Givens angles. Angles are transformed to
    have support (-inf, inf) */
void mcmclib_Givens_rotations(gsl_matrix* A, const gsl_vector* alpha);

/** \brief Givens angles and eigenvalues representation of a positive definite matrix

    'alpha' and 'sigma' vectors are reparametrized to have support (-inf, inf).
    @param M result
    @param alpha_sigma real-valued vector. First N x (N - 1) / 2 elmts. are
     Givens angles. Last N eltmts are eigenvalues.
*/
void mcmclib_Givens_representation(gsl_matrix* M, const gsl_vector* alpha_sigma);

/** \brief Givens angles and sing. values representation of a non-symm. matrix

    'alpha' and 'sigma' vectors are reparametrized to have support (-inf, inf).
    @param M result
    @param alpha12_sigma real-valued vector. First N x (N - 1) / 2 elmts. are
     Givens angles of matrix A, second N x (N - 1) / 2 elmts. are
     Givens angles of matrix B, Last N eltmts are singular values.
*/
void mcmclib_Givens_representation_asymm(gsl_matrix* M, const gsl_vector* alpha12_sigma);

/**@}*/
#endif
