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
#include <gsl/gsl_linalg.h>
#include "matrix.h"
#include "mcar_tilde.h"

mcmclib_mcar_tilde_lpdf* mcmclib_mcar_tilde_lpdf_alloc(size_t p, gsl_matrix* M) {
  mcmclib_mcar_tilde_lpdf* a = (mcmclib_mcar_tilde_lpdf*) malloc(sizeof(mcmclib_mcar_tilde_lpdf));
  assert(p>0);
  assert(M->size1 == M->size2);
  const size_t n = M->size1;
  const size_t offset = p * (p-1) / 2;
  a->p = p;
  a->n = n;
  a->B_tilde = gsl_matrix_alloc(p, p);
  a->Gamma = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(a->Gamma);
  a->alpha12sigma = gsl_vector_alloc(2*offset + p);
  gsl_vector_set_zero(a->alpha12sigma);
  a->M = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(a->M, M);
  a->m = gsl_vector_alloc(n);
  for(size_t i=0; i<n; i++) {
    size_t count = 0;
    for(size_t j=0; j<n; j++)
      count += gsl_matrix_get(M, i, j) == 1.0;
    gsl_vector_set(a->m, i, (double) count);
  }
  a->alphasigmag = gsl_vector_alloc(p * (p-1)/2 + p);
  gsl_vector_set_all(a->alphasigmag, -1.0);

  a->vcov = gsl_matrix_alloc(p*n, p*n);
  gsl_matrix_set_identity(a->vcov);

  a->Lambda_ij = gsl_matrix_alloc(p, p);
  a->Gammai = gsl_matrix_alloc(p, p);
  a->Block = gsl_matrix_alloc(p, p);

  mcmclib_mcar_tilde_lpdf_update_vcov(a);
  return a;
}

void mcmclib_mcar_tilde_lpdf_free(mcmclib_mcar_tilde_lpdf* p) {
  gsl_matrix_free(p->vcov);
  gsl_vector_free(p->m);
  gsl_matrix_free(p->M);
  gsl_vector_free(p->alpha12sigma);
  gsl_matrix_free(p->B_tilde);
  gsl_matrix_free(p->Gamma);
  gsl_vector_free(p->alphasigmag);

  gsl_matrix_free(p->Lambda_ij);
  gsl_matrix_free(p->Gammai);
  gsl_matrix_free(p->Block);
  free(p);
}

void mcmclib_mcar_tilde_lpdf_update_B_tilde(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_Givens_representation_asymm(p->B_tilde, p->alpha12sigma);
}

static int is_positive_definite(mcmclib_mcar_tilde_lpdf* p) {
  size_t offset = p->p;
  offset = offset * (offset-1) / 2;
  double sigma_0 = gsl_vector_get(p->alpha12sigma, 2*offset);
  if(sigma_0 > 0.0)
    return 0;
  return 1;
}

static void get_Lambda_L(gsl_matrix* Lambda_L, const gsl_matrix* A,
			 const gsl_matrix* B_tilde) {
  const size_t p = A->size1;
  gsl_matrix_set_zero(Lambda_L);
  gsl_matrix* AB = gsl_matrix_alloc(p, p);
  gsl_matrix_memcpy(AB, B_tilde);
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, A, AB);
  gsl_blas_dtrsm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, A, AB);
  gsl_matrix_memcpy(Lambda_L, AB);
  gsl_matrix_free(AB);
}

static inline void block_memcpy(gsl_matrix* dest, size_t i, size_t j, const gsl_matrix* src) {
  gsl_matrix_view block_view = gsl_matrix_submatrix(dest, i, j, src->size1, src->size2);
  gsl_matrix_memcpy(&block_view.matrix, src);
}

int mcmclib_mcar_tilde_lpdf_update_blocks(mcmclib_mcar_tilde_lpdf* p) {
  gsl_matrix* Lambda_ij = p->Lambda_ij;
  gsl_matrix* Gammai = p->Gammai;
  gsl_matrix* Block = p->Block;

  const size_t n = p->n;
  gsl_matrix* A = gsl_matrix_alloc(p->p, p->p);
  gsl_matrix_memcpy(A, p->Gamma);
  int status = mcmclib_cholesky_decomp(A);
  if(status) {
    gsl_matrix_free(A);
    return status;
  }
  gsl_matrix* A1 = gsl_matrix_alloc(p->p, p->p);
  gsl_matrix_memcpy(A1, A);
  gsl_linalg_cholesky_invert(A1);
  for(size_t i=0; i<(p->p - 1); i++)
    for (size_t j=i+1; j < p->p; j++)
      gsl_matrix_set(A, i, j, 0.0);

  get_Lambda_L(Lambda_ij, A, p->B_tilde);
  gsl_matrix_free(A);

  gsl_matrix_set_zero(p->vcov);
  for(size_t i=0; i<n; i++) {
    double mi = gsl_vector_get(p->m, i);
    gsl_matrix_memcpy(Gammai, A1);
    gsl_matrix_scale(Gammai, mi);
    for(size_t j=i; j<n; j++) {
      if(i == j) {
	block_memcpy(p->vcov, i * p->p, j * p->p, Gammai);
      } else if (gsl_matrix_get(p->M, i, j) == 1.0) {
	gsl_matrix_set_zero(Block);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0 / mi,
		       Gammai, Lambda_ij, 0.0, Block);
		       block_memcpy(p->vcov, j * p->p, i * p->p, Block);
      }
    }
  }

  gsl_matrix_free(A1);
  return GSL_SUCCESS;
}

int mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_mcar_tilde_lpdf_update_Gamma(p);
  mcmclib_mcar_tilde_lpdf_update_B_tilde(p);
  return mcmclib_mcar_tilde_lpdf_update_blocks(p);
}

double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mcar_tilde_lpdf* p = (mcmclib_mcar_tilde_lpdf*) in_p;
  if(!is_positive_definite(p))
    return log(0.0);
  int status = mcmclib_mcar_tilde_lpdf_update_vcov(p);
  if(status)
    return log(0.0);
  return mcmclib_mvnormzp_lpdf(p->vcov, x);
}

void mcmclib_mcar_tilde_lpdf_update_Gamma(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_Givens_representation(p->Gamma, p->alphasigmag);
}
