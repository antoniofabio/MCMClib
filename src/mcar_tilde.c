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
  int offset = p->p;
  offset = offset * (offset-1) / 2;
  double sigma_0 = gsl_vector_get(p->alpha12sigma, 2*offset);
  if(sigma_0 > 0.0)
    return 0;
  return 1;
}

static void get_Lambda_LU(gsl_matrix* Lambda_LU, int flag,
			  gsl_matrix* A, gsl_matrix* A1,
			  gsl_matrix* B_tilde) {
  int p = A->size1;
  gsl_matrix_set_zero(Lambda_LU);
  gsl_matrix* AB = gsl_matrix_alloc(p, p);
  gsl_matrix_set_zero(AB);
  if(flag) {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B_tilde, 0.0, AB);
  } else {
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, B_tilde, 0.0, AB);
  }
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AB, A1, 0.0, Lambda_LU);
  gsl_matrix_free(AB);
}

static inline void block_memcpy(gsl_matrix* dest, int i, int j, const gsl_matrix* src) {
  gsl_matrix_view block_view = gsl_matrix_submatrix(dest, i, j,
						    src->size1, src->size2);
  gsl_matrix_memcpy(&block_view.matrix, src);
}

int mcmclib_mcar_tilde_lpdf_update_blocks(mcmclib_mcar_tilde_lpdf* p) {
  gsl_matrix* Lambda_ij = p->Lambda_ij;
  gsl_matrix* Gammai = p->Gammai;
  gsl_matrix* Block = p->Block;

  int n = p->n;
  mcmclib_mcar_tilde_lpdf_update_B_tilde(p);
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
  for(int i=0; i<(p->p - 1); i++)
    for (int j=i+1; j < p->p; j++)
      gsl_matrix_set(A, i, j, 0.0);

  gsl_matrix* Lambda_L = gsl_matrix_alloc(p->p, p->p);
  get_Lambda_LU(Lambda_L, 0, A, A1, p->B_tilde);
  gsl_matrix* Lambda_U = gsl_matrix_alloc(p->p, p->p);
  get_Lambda_LU(Lambda_U, 1, A, A1, p->B_tilde);
  gsl_matrix_free(A);

  gsl_matrix_set_zero(p->vcov);
  for(int i=0; i<n; i++) {
    double mi = gsl_vector_get(p->m, i);
    gsl_matrix_memcpy(Gammai, A1);
    gsl_matrix_scale(Gammai, mi);
    for(int j=i; j<n; j++) {
      if(i == j) {
	block_memcpy(p->vcov, i * p->p, j * p->p, Gammai);
      } else if (gsl_matrix_get(p->M, i, j) == 1.0) {
	gsl_matrix_set_zero(Block);
	gsl_matrix_memcpy(Lambda_ij, Lambda_U);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0 / mi,
		       Gammai, Lambda_ij, 0.0, Block);
	block_memcpy(p->vcov, i * p->p, j * p->p, Block);
	gsl_matrix_set_zero(Block);
	gsl_matrix_memcpy(Lambda_ij, Lambda_L);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0 / mi,
		       Gammai, Lambda_ij, 0.0, Block);
	block_memcpy(p->vcov, j * p->p, i * p->p, Block);
      }
    }
  }

  gsl_matrix_free(A1);
  gsl_matrix_free(Lambda_L);
  gsl_matrix_free(Lambda_U);
  return GSL_SUCCESS;
}

int mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p) {
  int status = mcmclib_mcar_tilde_lpdf_update_blocks(p);
  if(status)
    return status;
  return mcmclib_cholesky_inverse(p->vcov);
}

/* Normal distrib. with 0 mean, no normalizing factors
   P is assumed to be the cholesky dec. of the precision matrix.
 */
static inline double _mvnorm(const gsl_matrix* P, const gsl_vector* e,
			     double ldet, gsl_vector* work) {
  gsl_vector_set_zero(work);
  gsl_blas_dsymv(CblasUpper, 1.0, P, e, 0.0, work);
  double ans = 0.0;
  gsl_blas_ddot(e, work, &ans);
  return -0.5 * ( ans + ldet );
}

double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mcar_tilde_lpdf* p = (mcmclib_mcar_tilde_lpdf*) in_p;
  if(!is_positive_definite(p))
    return log(0.0);
  mcmclib_mcar_tilde_lpdf_update_B_tilde(p);
  mcmclib_mcar_tilde_lpdf_update_Gamma(p);
  int status = mcmclib_mcar_tilde_lpdf_update_blocks(p);
  if(status)
    return log(0.0);
  gsl_vector* work = gsl_vector_alloc(x->size);
  gsl_matrix* tmp = gsl_matrix_alloc(x->size, x->size);
  gsl_matrix_memcpy(tmp, p->vcov);
  mcmclib_cholesky_decomp(tmp);
  double ldet = mcmclib_matrix_logtrace(tmp);
  gsl_matrix_free(tmp);
  double ans = _mvnorm(p->vcov, x, ldet, work);
  gsl_vector_free(work);
  return ans;
}

void mcmclib_mcar_tilde_lpdf_update_Gamma(mcmclib_mcar_tilde_lpdf* p) {
  mcmclib_Givens_representation(p->Gamma, p->alphasigmag);
}
