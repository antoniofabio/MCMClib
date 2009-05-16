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
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_blas.h>
#include "matrix.h"
#include "mvnorm.h"
#include "pois_model.h"

mcmclib_pois_model* mcmclib_pois_model_alloc(const gsl_matrix* X, const gsl_vector* y) {
  mcmclib_pois_model* a = (mcmclib_pois_model*) malloc(sizeof(mcmclib_pois_model));
  const int n = X->size1;
  const int p = X->size2;
  assert(y->size == n);
  a->X = gsl_matrix_alloc(n, p);
  gsl_matrix_memcpy(a->X, X);
  a->y = y;
  a->beta = gsl_vector_alloc(p);
  gsl_vector_set_all(a->beta, 0.0);
  a->b0 = gsl_vector_alloc(p);
  gsl_vector_set_zero(a->b0);
  a->B0 = gsl_matrix_alloc(p, p);
  a->mu = gsl_vector_alloc(n);
  a->work1 = gsl_vector_alloc(p);
  a->work2 = gsl_vector_alloc(p);
  a->offset = NULL;
  gsl_matrix* tmp = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(tmp);
  gsl_matrix_scale(tmp, 0.01);
  mcmclib_pois_model_set_prior_var(a, tmp);
  gsl_matrix_free(tmp);
  return a;
}

void mcmclib_pois_model_free(mcmclib_pois_model* p) {
  gsl_vector_free(p->work2);
  gsl_vector_free(p->work1);
  gsl_vector_free(p->mu);
  gsl_matrix_free(p->B0);
  gsl_vector_free(p->b0);
  gsl_vector_free(p->beta);
  gsl_matrix_free(p->X);
  free(p);
}

int mcmclib_pois_model_set_prior_mean(mcmclib_pois_model* p, const gsl_vector* b0) {
  gsl_vector_memcpy(p->b0, b0);
  return GSL_SUCCESS;
}
int mcmclib_pois_model_set_prior_var(mcmclib_pois_model* p, const gsl_matrix* B0) {
  gsl_matrix_memcpy(p->B0, B0);
  gsl_linalg_cholesky_decomp(p->B0);
  p->ldet = -mcmclib_matrix_logtrace(p->B0);
  gsl_matrix_memcpy(p->B0, B0);
  return GSL_SUCCESS;
}
int mcmclib_pois_model_set_offset(mcmclib_pois_model* p, const gsl_vector* offset) {
  assert(offset->size == p->X->size1);
  p->offset = offset;
  return GSL_SUCCESS;
}

static double lpois(double k, double lmu) {
  return lmu * k - gsl_sf_lnfact((int) k) - exp(lmu);
}

double mcmclib_pois_model_llik(mcmclib_pois_model* p, gsl_vector* x) {
  double ans = 0.0;
  gsl_vector_set_zero(p->mu);
  gsl_blas_dgemv(CblasNoTrans, 1.0, p->X, p->beta, 0.0, p->mu);
  if(p->offset)
    gsl_vector_add(p->mu, p->offset);
  for(int i=0; i < p->X->size1; i++)
    ans += lpois(gsl_vector_get(p->y, i), gsl_vector_get(p->mu, i));
  return ans;
}

double mcmclib_pois_model_lprior(mcmclib_pois_model* p, gsl_vector* x) {
  return mcmclib_mvnorm_lpdf_noinv(p->b0, p->B0, x, p->ldet, p->work1, p->work2);
}

double mcmclib_pois_model_lpdf(mcmclib_pois_model* p, gsl_vector* x) {
  return mcmclib_pois_model_lprior(p, x) + mcmclib_pois_model_llik(p, x);
}
