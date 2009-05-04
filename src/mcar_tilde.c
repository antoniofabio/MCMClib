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

mcmclib_mcar_tilde_lpdf* mcmclib_mcar_tilde_lpdf_alloc(int p, gsl_matrix* M) {
  mcmclib_mcar_tilde_lpdf* a = (mcmclib_mcar_tilde_lpdf*) malloc(sizeof(mcmclib_mcar_tilde_lpdf));
  assert(p>0);
  assert(M->size1 == M->size2);
  int n = M->size1;
  int offset = p * (p-1) / 2;
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
  for(int i=0; i<n; i++) {
    int count = 0;
    for(int j=0; j<n; j++)
      count += gsl_matrix_get(M, i, j) == 1.0;
    gsl_vector_set(a->m, i, (double) count);
  }
  a->alphasigmag = gsl_vector_alloc(p * (p-1)/2 + p);
  gsl_vector_set_all(a->alphasigmag, 0.0);

  a->vcov = gsl_matrix_alloc(p*n, p*n);
  gsl_matrix_set_identity(a->vcov);
  a->mu = gsl_vector_alloc(p*n);
  gsl_vector_set_zero(a->mu);
  a->mvnorm = mcmclib_mvnorm_lpdf_alloc(a->mu, a->vcov->data);

  a->Lambda_ij = gsl_matrix_alloc(p, p);
  a->Gammai = gsl_matrix_alloc(p, p);
  a->Block = gsl_matrix_alloc(p, p);

  mcmclib_mcar_tilde_lpdf_update_vcov(a);
  return a;
}

void mcmclib_mcar_tilde_lpdf_free(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_mvnorm_lpdf_free(p->mvnorm);
  gsl_matrix_free(p->vcov);
  gsl_vector_free(p->mu);
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
  int offset = p->p;
  offset = offset * (offset-1) / 2;
  double sigma_0 = gsl_vector_get(p->alpha12sigma, 2*offset);
  if((sigma_0 <= 0.0) || (sigma_0 >= 1.0))
    return 0;
  return 1;
}

static void get_inverse(gsl_matrix* A) {
  gsl_linalg_cholesky_decomp(A);
  gsl_linalg_cholesky_invert(A);
}

void mcmclib_matrix_inverse(gsl_matrix* A) {
  gsl_permutation* p = gsl_permutation_alloc(A->size1);
  gsl_matrix* A1 = gsl_matrix_alloc(A->size1, A->size1);
  int tmp=0;
  gsl_matrix_memcpy(A1, A);
  gsl_linalg_LU_decomp(A1, p, &tmp);
  gsl_linalg_LU_invert(A1, p, A);
  gsl_matrix_free(A1);
  gsl_permutation_free(p);
}

static void get_Lambda_LU(gsl_matrix* Lambda_LU, int flag, gsl_matrix* A, gsl_matrix* B_tilde) {
  int p = A->size1;
  gsl_matrix_set_zero(Lambda_LU);
  gsl_matrix* AB = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(AB);
  if(flag)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B_tilde, 0.0, AB);
  else
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, B_tilde, 0.0, AB);
  gsl_matrix* A1 = gsl_matrix_alloc(p, p);
  gsl_matrix_memcpy(A1, A);
  mcmclib_matrix_inverse(A1);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AB, A1, 0.0, Lambda_LU);
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
  mcmclib_mcar_tilde_lpdf_update_B_tilde(p);
  gsl_matrix* A = gsl_matrix_alloc(p->p, p->p);
  gsl_matrix_memcpy(A, p->Gamma);
  gsl_linalg_cholesky_decomp(A);
  for(int i=0; i<(p->p - 1); i++)
    for (int j=i+1; j < p->p; j++)
      gsl_matrix_set(A, i, j, 0.0);

  gsl_matrix* Lambda_L = gsl_matrix_alloc(p->p, p->p);
  get_Lambda_LU(Lambda_L, 0, A, p->B_tilde);
  gsl_matrix* Lambda_U = gsl_matrix_alloc(p->p, p->p);
  get_Lambda_LU(Lambda_U, 1, A, p->B_tilde);

  for(int i=0; i<n; i++) {
    gsl_matrix_memcpy(Gammai, p->Gamma);
    get_inverse(Gammai);
    gsl_matrix_scale(Gammai, gsl_vector_get(p->m, i));
    for(int j=0; j<n; j++) {
      gsl_matrix_set_zero(Block);
      if((gsl_matrix_get(p->M, i, j) == 1.0) || (i==j)) {
	if(i<j) {
	  gsl_matrix_memcpy(Lambda_ij, Lambda_U);
	  gsl_matrix_scale(Lambda_ij, 1.0 / gsl_vector_get(p->m, i));
	} else if(i>j) {
	  gsl_matrix_memcpy(Lambda_ij, Lambda_L);
	  gsl_matrix_scale(Lambda_ij, 1.0 / gsl_vector_get(p->m, j));
	} else {
	  gsl_matrix_set_identity(Lambda_ij);
	  gsl_matrix_scale(Lambda_ij, -1.0);
	}
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, Gammai, Lambda_ij, 0.0, Block);
      }
      block_memcpy(p->vcov, i * p->p, j * p->p, Block);
    }
  }

  gsl_matrix_free(Lambda_L);
  gsl_matrix_free(Lambda_U);
  gsl_matrix_free(A);
}

void mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_mcar_tilde_lpdf_update_blocks(p);
  mcmclib_matrix_inverse(p->vcov);
}

double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mcar_tilde_lpdf* p = (mcmclib_mcar_tilde_lpdf*) in_p;
  if(!is_positive_definite(p))
    return log(0.0);
  mcmclib_mcar_tilde_lpdf_update_B_tilde(p);
  mcmclib_mcar_tilde_lpdf_update_Gamma(p);
  mcmclib_mcar_tilde_lpdf_update_vcov(p);
  return mcmclib_mvnorm_lpdf_compute(p->mvnorm, x);
}

void mcmclib_mcar_tilde_lpdf_update_Gamma(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_Givens_representation(p->Gamma, p->alphasigmag);
}
