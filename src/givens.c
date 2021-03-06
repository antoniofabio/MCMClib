/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include "givens.h"

/*sum of the numbers from 'b' to 'a'*/
static size_t partsum(size_t a, size_t b) {
  if(b>a)
    return 0;
  return (a*(a+1) - b*(b-1))/2;
}

static size_t ALPHA(size_t i, size_t j, size_t p) {
  assert(i < j);
  return partsum(p-1, p-i) + j - i - 1;
}

static double alpha_get(const gsl_vector* a, size_t i, size_t j, size_t p) {
  return gsl_vector_get(a, ALPHA(i, j, p));
}

static void Givens_set_Shij(gsl_matrix* S, size_t i, size_t j, double alpha_ij) {
  gsl_matrix_set_identity(S);
  gsl_matrix_set(S, i, i, cos(alpha_ij));
  gsl_matrix_set(S, j, j, cos(alpha_ij));
  gsl_matrix_set(S, i, j, sin(alpha_ij));
  gsl_matrix_set(S, j, i, -sin(alpha_ij));
}

void mcmclib_Givens_rotations(gsl_matrix* A, const gsl_vector* alpha) {
  assert(A->size1 == A->size2);
  const size_t p = A->size1;
  assert(alpha->size == (p * (p-1) / 2));
  gsl_matrix* S[2];
  for(int h=0; h<2; h++)
    S[h] = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(S[1]);

  for(size_t i=0; i<(p-1); i++) {
    for(size_t j=i+1; j<p; j++) {
      double bi = exp(alpha_get(alpha, i, j, p));
      bi = M_PI_2 * (bi - 1.0) / (bi + 1.0);
      Givens_set_Shij(S[0], i, j, bi);
      gsl_matrix_set_zero(A);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S[1], S[0], 0.0, A);
      gsl_matrix_memcpy(S[1], A);
    }
  }
  gsl_matrix_memcpy(A, S[1]);

  for(int h=0; h<2; h++)
    gsl_matrix_free(S[h]);
}

static void vSortDesc(gsl_vector* v) {
  gsl_vector_scale(v, -1.0);
  gsl_sort_vector(v);
  gsl_vector_scale(v, -1.0);
}

void mcmclib_Givens_representation(gsl_matrix* M,
				   const gsl_vector* alpha_sigma) {
  size_t n = M->size1;
  size_t offset = n * (n-1) / 2;

  gsl_vector* alpha1 = gsl_vector_alloc(offset);
  for(size_t i=0; i<offset; i++)
    gsl_vector_set(alpha1, i, gsl_vector_get(alpha_sigma, i));
  gsl_vector* sigma1 = gsl_vector_alloc(n);
  for(size_t i=0; i<n; i++) {
    double bi = exp(gsl_vector_get(alpha_sigma, i + offset));
    gsl_vector_set(sigma1, i, bi);
  }
  vSortDesc(sigma1);

  gsl_matrix* A = gsl_matrix_alloc(n, n);
  mcmclib_Givens_rotations(A, alpha1);
  gsl_matrix* D = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(D);
  for(size_t i=0; i<n; i++)
    gsl_matrix_set(D, i, i, gsl_vector_get(sigma1, i));
  gsl_matrix* AD = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(AD);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, D, 0.0, AD);
  gsl_matrix_set_zero(M);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, AD, A, 0.0, M);

  gsl_matrix_free(AD);
  gsl_matrix_free(D);
  gsl_matrix_free(A);
  gsl_vector_free(alpha1);
  gsl_vector_free(sigma1);
}

/* rebuild asymm. matrix from its SVD decomposition */
static void anti_SVD(gsl_matrix* A,
		     const gsl_matrix* P1,
		     const gsl_matrix* P2,
		     const gsl_vector* sigma) {
  const size_t p = A->size1;
  gsl_matrix* Delta = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(Delta);
  for(size_t i=0; i<p; i++)
    gsl_matrix_set(Delta, i, i, exp(gsl_vector_get(sigma, i)));
  gsl_matrix* P1Delta = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(P1Delta);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P1, Delta, 0.0, P1Delta);
  gsl_matrix_set_zero(A);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P1Delta, P2, 0.0, A);
  gsl_matrix_free(Delta);
  gsl_matrix_free(P1Delta);
}

void mcmclib_Givens_representation_asymm(gsl_matrix* M, const gsl_vector* alpha12_sigma) {
  const size_t n = M->size1;
  const size_t offset = n*(n-1)/2;
  gsl_matrix* P1 = gsl_matrix_alloc(n, n);
  gsl_matrix* P2 = gsl_matrix_alloc(n, n);
  gsl_vector_const_view alpha1 = gsl_vector_const_subvector(alpha12_sigma, 0, offset);
  gsl_vector_const_view alpha2 = gsl_vector_const_subvector(alpha12_sigma, offset, offset);
  gsl_vector_const_view sigma = gsl_vector_const_subvector(alpha12_sigma, 2*offset, n);
  mcmclib_Givens_rotations(P1, &alpha1.vector);
  mcmclib_Givens_rotations(P2, &alpha2.vector);
  gsl_vector* sigmas = gsl_vector_alloc(n);
  gsl_vector_memcpy(sigmas, &sigma.vector);
  vSortDesc(sigmas);
  anti_SVD(M, P1, P2, sigmas);
  gsl_vector_free(sigmas);
  gsl_matrix_free(P1);
  gsl_matrix_free(P2);
}
