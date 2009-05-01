/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <gsl/gsl_sort_vector.h>
#include "mcar_model.h"

mcmclib_mcar_model* mcmclib_mcar_model_alloc(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e) {
  mcmclib_mcar_model* a = (mcmclib_mcar_model*) malloc(sizeof(mcmclib_mcar_model));
  a->lpdf = m;
  a->e = e;

  int p = m->p;
  gsl_matrix* V = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(V);
  a->w = mcmclib_wishart_lpdf_alloc(V, p);
  gsl_matrix_free(V);
  return a;
}

void mcmclib_mcar_model_free(mcmclib_mcar_model* p) {
  mcmclib_wishart_lpdf_free(p->w);
  free(p);
}

double mcmclib_mcar_model_alpha1_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha1) {
  gsl_vector* tmp = gsl_vector_alloc(alpha1->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha1);
  gsl_vector* alpha1b = gsl_vector_alloc(alpha1->size);
  for(int i=0; i < alpha1->size; i++) {
    double bi = exp(gsl_vector_get(alpha1, i));
    gsl_vector_set(alpha1b, i, M_PI_2 * (bi - 1.0) / (bi + 1.0));
  }
  gsl_vector_memcpy(p->lpdf->alpha1, alpha1b);
  gsl_vector_free(alpha1b);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha1, tmp);
  gsl_vector_free(tmp);
  return ans;
}

double mcmclib_mcar_model_alpha2_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha2) {
  gsl_vector* tmp = gsl_vector_alloc(alpha2->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha2);
  gsl_vector* alpha2b = gsl_vector_alloc(alpha2->size);
  for(int i=0; i < alpha2->size; i++) {
    double bi = exp(gsl_vector_get(alpha2, i));
    gsl_vector_set(alpha2b, i, M_PI_2 * (bi - 1.0) / (bi + 1.0));
  }
  gsl_vector_memcpy(p->lpdf->alpha2, alpha2b);
  gsl_vector_free(alpha2b);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha2, tmp);
  gsl_vector_free(tmp);
  return ans;
}

double mcmclib_mcar_model_sigma_lpdf(mcmclib_mcar_model* p, gsl_vector* sigma) {
  int P = sigma->size;
  gsl_vector* tmp = gsl_vector_alloc(P);
  gsl_vector_memcpy(tmp, p->lpdf->sigma);
  gsl_vector* sigma1 = gsl_vector_alloc(P);
  gsl_vector_memcpy(sigma1, sigma);
  gsl_vector_scale(sigma1, -1.0);
  gsl_sort_vector(sigma1);
  gsl_vector_scale(sigma1, -1.0);
  for(int i=0; i<P; i++) {
    double bs = exp(gsl_vector_get(sigma1, i));
    gsl_vector_set(sigma1, i, 0.5 * (bs - 1.0) / (bs + 1.0) + 0.5 );
  }
  gsl_vector_memcpy(p->lpdf->sigma, sigma1);
  gsl_vector_free(sigma1);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->sigma, tmp);
  gsl_vector_free(tmp);
  return ans;
}

/* Givens angles and eigenvalues representation of a pos.def. matrix.
   Result goes in M */
void givens_representation(gsl_matrix* M, gsl_vector* alpha, gsl_vector* sigma) {
  int n = M->size1;

  gsl_matrix* A = gsl_matrix_alloc(n, n);
  mcmclib_Givens_rotations(A, alpha);
  gsl_matrix* D = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(D);
  for(int i=0; i<n; i++)
    gsl_matrix_set(D, i, i, gsl_vector_get(sigma, i));
  gsl_matrix* AD = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(AD);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, D, 0.0, AD);
  gsl_matrix_set_zero(M);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, AD, A, 0.0, M);

  gsl_matrix_free(AD);
  gsl_matrix_free(D);
  gsl_matrix_free(A);
}

double mcmclib_mcar_model_alphasigma_lpdf(mcmclib_mcar_model* p, gsl_vector* alphasigma) {
  int P = p->lpdf->p;
  gsl_matrix* Gamma = gsl_matrix_alloc(P, P);
  gsl_vector_view alpha_v = gsl_vector_subvector(alphasigma, 0, P*(P-1)/2);
  gsl_vector* alpha = &alpha_v.vector;
  gsl_vector* alpha1 = gsl_vector_alloc(P*(P-1)/2);
  for(int i=0; i<P*(P-1)/2; i++) {
    double bi = exp(gsl_vector_get(alpha, i));
    gsl_vector_set(alpha1, i, M_PI_2 * (bi - 1.0) / (bi + 1.0));
  }
  gsl_vector_view sigma_v = gsl_vector_subvector(alphasigma, P*(P-1)/2, P);
  gsl_vector* sigma = &sigma_v.vector;
  gsl_vector* sigma1 = gsl_vector_alloc(P);
  for(int i=0; i<P; i++)
    gsl_vector_set(sigma1, i, exp(gsl_vector_get(sigma, i)));
  givens_representation(Gamma, alpha, sigma1);
  gsl_vector_free(sigma1);
  gsl_vector_view gamma_v = gsl_vector_view_array(Gamma->data, P*(P-1)/2);
  double ans = mcmclib_mcar_model_Gamma_lpdf(p, &gamma_v.vector);
  gsl_matrix_free(Gamma);
  return ans;
}

double mcmclib_mcar_model_Gamma_lpdf(mcmclib_mcar_model* p, gsl_vector* gamma) {
  int P = p->lpdf->p;
  gsl_matrix* tmp = gsl_matrix_alloc(P, P);
  gsl_matrix_view gamma_v = gsl_matrix_view_vector(gamma, P, P);
  gsl_matrix_memcpy(tmp, p->lpdf->Gamma);
  gsl_matrix_memcpy(p->lpdf->Gamma, &(gamma_v.matrix));
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e)
    + mcmclib_wishart_lpdf_compute(p->w, gamma);
  gsl_matrix_memcpy(p->lpdf->Gamma, tmp);
  gsl_matrix_free(tmp);
  return ans;
}
