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
#include <gsl/gsl_randist.h>
#include "matrix.h"
#include "mcar_model.h"

/* tolerance for pos. def. conditions */
#define TOL 5e-2

mcmclib_mcar_model* mcmclib_mcar_model_alloc(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e) {
  mcmclib_mcar_model* a = (mcmclib_mcar_model*) malloc(sizeof(mcmclib_mcar_model));
  a->lpdf = m;
  a->e = e;

  int p = m->p;
  gsl_matrix* V = gsl_matrix_alloc(p, p);
  gsl_matrix_set_identity(V);
  a->w = mcmclib_iwishart_lpdf_alloc(V, p);
  gsl_matrix_free(V);
  return a;
}

void mcmclib_mcar_model_free(mcmclib_mcar_model* p) {
  mcmclib_iwishart_lpdf_free(p->w);
  free(p);
}

static double alphaij_logderiv(double x) {
  return x - 2 * log(exp(x) + 1.0);
}

static double alpha12sigma_logderiv(int p, const gsl_vector* x) {
  double ans = 0.0;
  const int offset = p*(p-1);
  for(int n=0; n < offset; n++) {
    double xn = gsl_vector_get(x, n);
    if(abs(xn) >= log(M_PI / TOL - 1.0))
      return log(0.0);
    ans += alphaij_logderiv(xn);
  }
  for(int n=offset; n < x->size; n++) {
    double xn = gsl_vector_get(x, n);
    if((xn <= log(TOL)) || (xn >= log(1.0 - TOL)))
      return log(0.0);
    ans += xn;
  }
  return ans;
}

double mcmclib_mcar_model_alpha12sigma_lpdf(void* in_p, gsl_vector* alpha12sigma) {
  mcmclib_mcar_model* p = (mcmclib_mcar_model*) in_p;
  int n = p->lpdf->p;
  gsl_vector_view s_v = gsl_vector_subvector(alpha12sigma, n*(n-1), n);
  if(!mcmclib_vector_is_sorted_desc(&s_v.vector))
    return log(0.0);
  double prior = alpha12sigma_logderiv(n, alpha12sigma);
  if(!gsl_finite(prior))
    return prior;
  gsl_vector* tmp = gsl_vector_alloc(alpha12sigma->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha12sigma);
  gsl_vector_memcpy(p->lpdf->alpha12sigma, alpha12sigma);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha12sigma, tmp);
  gsl_vector_free(tmp);
  return ans + prior;
}

static double alphasigma_logderiv(int p, const gsl_vector* x) {
  double ans = 0.0;
  int offset = p*(p-1)/2;
  for(int n=0; n < offset; n++) {
    double xn = gsl_vector_get(x, n);
    if(abs(xn) >= log(M_PI / TOL - 1.0))
      return log(0.0);
    ans += alphaij_logderiv(xn);
  }
  for(int n=offset; n < x->size; n++) {
    double xn = gsl_vector_get(x, n);
    if(xn <= log(TOL))
      return log(0.0);
    ans += 0;
  }
  return ans;
}

double mcmclib_mcar_model_alphasigma_lpdf(void* in_p, gsl_vector* alphasigma) {
  mcmclib_mcar_model* p = (mcmclib_mcar_model*) in_p;

  const int n = p->lpdf->p;
  gsl_vector_view s_v = gsl_vector_subvector(alphasigma, n*(n-1)/2, n);
  if(!mcmclib_vector_is_sorted_desc(&s_v.vector))
    return log(0.0);

  double ans = alphasigma_logderiv(p->lpdf->p, alphasigma);
  if(!gsl_finite(ans))
    return ans;

  gsl_vector* tmp = gsl_vector_alloc(alphasigma->size);
  gsl_vector_memcpy(tmp, p->lpdf->alphasigmag);
  gsl_vector_memcpy(p->lpdf->alphasigmag, alphasigma);
  double lik = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  ans += lik;

  gsl_vector_memcpy(p->lpdf->alphasigmag, tmp);
  gsl_vector_free(tmp);
  return ans;
}

double mcmclib_mcar_model_phi_fcond(mcmclib_mcar_model* in_p, int i, gsl_vector* x) {
  mcmclib_mcar_tilde_lpdf* lpdf = in_p->lpdf;
  const int p = lpdf->p;
  if(x->size != p) {
    static char msg[1024];
    sprintf(msg, "'x' vector size is %ld, it should be %d", x->size, p);
    GSL_ERROR(msg, GSL_FAILURE);
  }
  assert(x->size == p);
  gsl_matrix* W = lpdf->M;
  const double mi = 1.0 / gsl_vector_get(lpdf->m, i);
  gsl_matrix* Lambdaij = lpdf->Lambda_ij;
  gsl_vector* mean = gsl_vector_alloc(p);
  gsl_vector_set_zero(mean);
  for(int j=0; j < lpdf->n; j++) {
    if(gsl_matrix_get(W, i, j) == 1.0) {
      gsl_vector_view phij_v = gsl_vector_subvector(in_p->e, j*p, p);
      gsl_blas_dgemv(CblasNoTrans, mi, Lambdaij, &phij_v.vector, 1.0, mean);
    }
  }
  gsl_matrix* Gammai = lpdf->Gammai;
  gsl_matrix_memcpy(Gammai, lpdf->Gamma);
  gsl_matrix_scale(Gammai, mi);
  mcmclib_mvnorm_lpdf* tmp = mcmclib_mvnorm_lpdf_alloc(mean, Gammai->data);
  double ans = mcmclib_mvnorm_lpdf_compute(tmp, x);
  gsl_vector_free(mean);
  mcmclib_mvnorm_lpdf_free(tmp);
  return ans;
}
