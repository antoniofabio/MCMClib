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

static int ALPHA(int i, int j, int p) {
  assert(i < j);
  if(i == 0)
    return (j - i) - 1;
  return p - i - 1 + ALPHA(i-1, j, p);
}

static double alpha_get(const gsl_vector* a, int i, int j) {
  return gsl_vector_get(a, ALPHA(i, j, a->size));
}

static void Givens_set_Shij(gsl_matrix* S, int i, int j, double alpha_ij) {
  gsl_matrix_set_identity(S);
  gsl_matrix_set(S, i, i, cos(alpha_ij));
  gsl_matrix_set(S, j, j, cos(alpha_ij));
  gsl_matrix_set(S, i, j, sin(alpha_ij));
  gsl_matrix_set(S, j, i, -sin(alpha_ij));
}

void mcmclib_Givens_rotations(gsl_matrix* A, const gsl_vector* alpha) {
  assert(A->size1 == A->size2);
  int p = A->size1;
  assert(alpha->size == (p * (p-1) / 2));
  gsl_matrix* S[2];
  for(int h=0; h<2; h++)
    S[h] = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(S[1]);

  for(int i=0; i<(p-1); i++) {
    for(int j=i+1; j<p; j++) {
      double bi = exp(alpha_get(alpha, i, j));
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
  int n = M->size1;
  int offset = n * (n-1) / 2;

  gsl_vector* alpha1 = gsl_vector_alloc(offset);
  for(int i=0; i<offset; i++)
    gsl_vector_set(alpha1, i, gsl_vector_get(alpha_sigma, i));
  gsl_vector* sigma1 = gsl_vector_alloc(n);
  for(int i=0; i<n; i++) {
    double bi = exp(gsl_vector_get(alpha_sigma, i + offset));
    gsl_vector_set(sigma1, i, bi);
  }
  vSortDesc(sigma1);

  gsl_matrix* A = gsl_matrix_alloc(n, n);
  mcmclib_Givens_rotations(A, alpha1);
  gsl_matrix* D = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(D);
  for(int i=0; i<n; i++)
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
static void anti_SVD(gsl_matrix* A, gsl_matrix* P1, gsl_matrix* P2, const gsl_vector* sigma) {
  int p = A->size1;
  gsl_matrix* Delta = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(Delta);
  for(int i=0; i<p; i++)
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
  int n = M->size1;
  int offset = n*(n-1)/2;
  gsl_matrix* P1 = gsl_matrix_alloc(n, n);
  gsl_matrix* P2 = gsl_matrix_alloc(n, n);
  gsl_vector_const_view alpha1 = gsl_vector_const_subvector(alpha12_sigma, 0, offset);
  gsl_vector_const_view alpha2 = gsl_vector_const_subvector(alpha12_sigma, offset, offset);
  gsl_vector_const_view sigma = gsl_vector_const_subvector(alpha12_sigma, 2*offset, n);
  gsl_vector* sigmas = gsl_vector_alloc(n);
  gsl_vector_memcpy(sigmas, &sigma.vector);
  vSortDesc(sigmas);
  mcmclib_Givens_rotations(P1, &alpha1.vector);
  mcmclib_Givens_rotations(P2, &alpha2.vector);
  anti_SVD(M, P1, P2, sigmas);
  gsl_vector_free(sigmas);
  gsl_matrix_free(P1);
  gsl_matrix_free(P2);
}
