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
#include <gsl/gsl_linalg.h>
#include "mcar_tilde.h"

mcmclib_mcar_tilde_lpdf* mcmclib_mcar_tilde_lpdf_alloc(int p, int n, gsl_matrix* M) {
  mcmclib_mcar_tilde_lpdf* a = (mcmclib_mcar_tilde_lpdf*) malloc(sizeof(mcmclib_mcar_tilde_lpdf));
  assert(p>0);
  assert(n>0);
  a->p = p;
  a->n = n;
  a->B_tilde = gsl_matrix_alloc(p, p);
  a->Gamma = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(a->Gamma);
  a->alpha1 = gsl_vector_alloc(p * (p-1) / 2);
  a->alpha2 = gsl_vector_alloc(p * (p-1) / 2);
  a->sigma = gsl_vector_alloc(p);
  gsl_vector_set_all(a->alpha1, 0.0);
  gsl_vector_set_all(a->alpha2, 0.0);
  gsl_vector_set_all(a->sigma, 0.4);
  gsl_vector_set(a->sigma, 0, 0.5);
  a->M = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(a->M, M);
  a->m = gsl_vector_alloc(n);
  for(int i=0; i<n; i++) {
    int count = 0;
    for(int j=0; j<n; j++)
      count += gsl_matrix_get(M, i, j) == 1.0;
    gsl_vector_set(a->m, i, (double) count);
  }

  a->vcov = gsl_matrix_alloc(p*n, p*n);
  gsl_matrix_set_identity(a->vcov);
  a->mu = gsl_vector_alloc(p*n);
  gsl_vector_set_zero(a->mu);
  a->mvnorm = mcmclib_mvnorm_lpdf_alloc(a->mu, a->vcov->data);

  a->Lambda_ij = gsl_matrix_alloc(p, p);
  a->Gammai = gsl_matrix_alloc(p, p);
  a->Block = gsl_matrix_alloc(p, p);

  return a;
}

void mcmclib_mcar_tilde_lpdf_free(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_mvnorm_lpdf_free(p->mvnorm);
  gsl_matrix_free(p->vcov);
  gsl_vector_free(p->mu);
  gsl_vector_free(p->m);
  gsl_matrix_free(p->M);
  gsl_vector_free(p->sigma);
  gsl_vector_free(p->alpha1);
  gsl_vector_free(p->alpha2);
  gsl_matrix_free(p->B_tilde);
  gsl_matrix_free(p->Gamma);

  gsl_matrix_free(p->Lambda_ij);
  gsl_matrix_free(p->Gammai);
  gsl_matrix_free(p->Block);
  free(p);
}

void mcmclib_mcar_tilde_lpdf_set_alpha(mcmclib_mcar_tilde_lpdf* p,
				       gsl_vector* alpha1, gsl_vector* alpha2) {
  gsl_vector_memcpy(p->alpha1, alpha1);
  gsl_vector_memcpy(p->alpha2, alpha2);
}

void mcmclib_mcar_tilde_lpdf_set_sigma(mcmclib_mcar_tilde_lpdf* p,
				       gsl_vector* sigma) {
  gsl_vector_memcpy(p->sigma, sigma);
}

static int ALPHA(int i, int j, int p) {
  assert(i < j);
  if(i == 0)
    return (j - i) - 1;
  return p - i - 1 + ALPHA(i-1, j, p);
}

static double alpha_get(gsl_vector* a, int i, int j) {
  return gsl_vector_get(a, ALPHA(i, j, a->size));
}

static void Givens_set_Shij(gsl_matrix* S, int i, int j, double alpha_ij) {
  gsl_matrix_set_identity(S);
  gsl_matrix_set(S, i, i, cos(alpha_ij));
  gsl_matrix_set(S, j, j, cos(alpha_ij));
  gsl_matrix_set(S, i, j, sin(alpha_ij));
  gsl_matrix_set(S, j, i, -sin(alpha_ij));
}

void mcmclib_Givens_rotations(gsl_matrix* A, gsl_vector* alpha) {
  assert(A->size1 == A->size2);
  int p = A->size1;
  assert(alpha->size == (p * (p-1) / 2));
  gsl_matrix* S[2];
  for(int h=0; h<2; h++)
    S[h] = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(A);
  gsl_matrix_set_identity(S[1]);

  for(int i=0; i<(p-1); i++) {
    for(int j=i+1; j<p; j++) {
      Givens_set_Shij(S[0], i, j, alpha_get(alpha, i, j));
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S[1], S[0], 0.0, A);
      gsl_matrix_memcpy(S[1], A);
    }
  }

  for(int h=0; h<2; h++)
    gsl_matrix_free(S[h]);
}

/* rebuild asymm. matrix from its SVD decomposition */
static void anti_SVD(gsl_matrix* A, gsl_matrix* P1, gsl_matrix* P2, gsl_vector* sigma) {
  int p = A->size1;
  gsl_matrix* Delta = gsl_matrix_alloc(p, p);
  gsl_matrix* P1Delta = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(Delta);
  for(int i=0; i<p; i++)
    gsl_matrix_set(Delta, i, i, gsl_vector_get(sigma, i));
  gsl_matrix_set_zero(P1Delta);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P1, Delta, 0.0, P1Delta);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P1Delta, P2, 0.0, A);
  gsl_matrix_free(Delta);
  gsl_matrix_free(P1Delta);
}

static void get_B_tilde(gsl_matrix* A, gsl_vector* sigma,
			gsl_vector* alpha1, gsl_vector* alpha2) {
  int p = A->size1;
  gsl_matrix* P1 = gsl_matrix_alloc(p, p);
  gsl_matrix* P2 = gsl_matrix_alloc(p, p);
  mcmclib_Givens_rotations(P1, alpha1);
  mcmclib_Givens_rotations(P2, alpha2);
  anti_SVD(A, P1, P2, sigma);
  gsl_matrix_free(P1);
  gsl_matrix_free(P2);
}

void mcmclib_mcar_tilde_lpdf_update_B_tilde(mcmclib_mcar_tilde_lpdf* p) {
  get_B_tilde(p->B_tilde, p->sigma, p->alpha1, p->alpha2);
}

static int is_positive_definite(mcmclib_mcar_tilde_lpdf* p) {
  double sigma_0 = gsl_vector_get(p->sigma, 0);
  if((sigma_0 <= 0.0) ||(sigma_0 >= 1.0))
    return 0;
  return 1;
}

static void get_inverse(gsl_matrix* A) {
  gsl_linalg_cholesky_decomp(A);
  gsl_linalg_cholesky_invert(A);
}

static void get_inverse_LU(gsl_matrix* A) {
  gsl_permutation* p = gsl_permutation_alloc(A->size1);
  gsl_matrix* A1 = gsl_matrix_alloc(A->size1, A->size1);
  int tmp=0;
  gsl_matrix_memcpy(A1, A);
  gsl_linalg_LU_decomp(A1, p, &tmp);
  gsl_linalg_LU_invert(A1, p, A);
  gsl_matrix_free(A1);
  gsl_permutation_free(p);
}

static void get_Lambda_ij(gsl_matrix* Lambda_ij, int i, int j,
			  gsl_vector* m, gsl_matrix* M,
			  gsl_matrix* A, gsl_matrix* B_tilde) {
  int p = A->size1;
  if(gsl_matrix_get(M, i, j) != 1.0) {
    gsl_matrix_set_zero(Lambda_ij);
    return;
  }
  gsl_matrix* AB = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(AB);
  if(i < j)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B_tilde, 0.0, AB);
  else
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, B_tilde, 0.0, AB);
  gsl_matrix* A1 = gsl_matrix_alloc(p, p);
  gsl_matrix_memcpy(A1, A);
  get_inverse(A1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, AB, A1, 0.0, Lambda_ij);
  gsl_matrix_scale(Lambda_ij, 1.0 / (i<j) ? gsl_vector_get(m, i) : gsl_vector_get(m, j) );
  gsl_matrix_free(AB);
  gsl_matrix_free(A1);
}

static void block_memcpy(gsl_matrix* dest, int i, int j, gsl_matrix* src) {
  int p = src->size1;
  int q = src->size2;
  for(int i1 = 0; i1 < p; i1++)
    for(int j1 = 0; j1 < q; j1++)
      gsl_matrix_set(dest, i1+i, j1+j, gsl_matrix_get(src, i1, j1));
}

void mcmclib_mcar_tilde_lpdf_update_blocks(mcmclib_mcar_tilde_lpdf* p) {
  gsl_matrix* Lambda_ij = p->Lambda_ij;
  gsl_matrix* Gammai = p->Gammai;
  gsl_matrix* Block = p->Block;

  int n = p->n;
  gsl_matrix* A = gsl_matrix_alloc(p->p, p->p);
  gsl_matrix_memcpy(A, p->Gamma);
  gsl_linalg_cholesky_decomp(A);
  for(int i=0; i<(p->p - 1); i++)
    for (int j=i+1; j < p->p; j++)
      gsl_matrix_set(A, i, j, gsl_matrix_get(A, j, i));

  gsl_matrix_set_zero(Block);
  for(int i=0; i<n; i++) {
    gsl_matrix_memcpy(Gammai, p->Gamma);
    gsl_matrix_scale(Gammai, 1.0 / gsl_vector_get(p->m, i));
    get_inverse(Gammai);
    for(int j=0; j<n; j++) {
      get_Lambda_ij(Lambda_ij, i, j, p->m, p->M, A, p->B_tilde);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Gammai, Lambda_ij, 0.0, Block);
      gsl_matrix_scale(Block, -1.0);
      block_memcpy(p->vcov, i * p->p, j * p->p, Block);
    }
  }

  gsl_matrix_free(A);
}

void mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_mcar_tilde_lpdf_update_blocks(p);
  get_inverse_LU(p->vcov);
}

double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mcar_tilde_lpdf* p = (mcmclib_mcar_tilde_lpdf*) in_p;
  if(!is_positive_definite(p))
    return log(0.0);
  mcmclib_mcar_tilde_lpdf_update_B_tilde(p);
  mcmclib_mcar_tilde_lpdf_update_vcov(p);
  return mcmclib_mvnorm_lpdf_compute(p->mvnorm, x);
}
