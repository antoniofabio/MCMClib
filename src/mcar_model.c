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

static double mvgauss(gsl_vector* x, double sigma) {
  double ans = 0.0;
  for(int i=0; i<x->size; i++)
    ans += log(gsl_ran_gaussian_pdf(gsl_vector_get(x, i), sigma));
  return ans;
}

static double alphaij_logderiv(double x) {
  return x - 2 * log(exp(x) + 1.0);
}

static double alpha12sigma_logderiv(int p, const gsl_vector* x) {
  double ans = 0.0;
  const int offset = p*(p-1);
  for(int n=offset; n < x->size; n++) {
    double xn = gsl_vector_get(x, n);
    if(xn >= 0.0)
      return log(0.0);
    ans += xn;
  }
  for(int n=0; n < offset; n++)
    ans += alphaij_logderiv(gsl_vector_get(x, n));
  return ans;
}

double mcmclib_mcar_model_alpha12sigma_lpdf(void* in_p, gsl_vector* alpha12sigma) {
  mcmclib_mcar_model* p = (mcmclib_mcar_model*) in_p;
  int n = p->lpdf->p;
  gsl_vector_view s_v = gsl_vector_subvector(alpha12sigma, n*(n-1), n);
  if(!mcmclib_vector_is_sorted_desc(&s_v.vector))
    return log(0.0);
  gsl_vector* tmp = gsl_vector_alloc(alpha12sigma->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha12sigma);
  gsl_vector_memcpy(p->lpdf->alpha12sigma, alpha12sigma);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha12sigma, tmp);
  gsl_vector_free(tmp);
  return ans + alpha12sigma_logderiv(n, alpha12sigma);
}

static double iwishart_alphasigma(mcmclib_iwishart_lpdf* p, gsl_vector* as) {
  int n = p->Psi->size1;
  gsl_vector_view s_v = gsl_vector_subvector(as, n*(n-1)/2, n);
  if(!mcmclib_vector_is_sorted_desc(&s_v.vector))
    return log(0.0);
  gsl_matrix* Gamma = gsl_matrix_alloc(n, n);
  mcmclib_Givens_representation(Gamma, as);
  gsl_vector_view gamma_v = gsl_vector_view_array(Gamma->data, n*n);
  double ans = mcmclib_iwishart_lpdf_compute(p, &gamma_v.vector);
  gsl_matrix_free(Gamma);
  return ans;
}

static double alphasigma_logderiv(int p, const gsl_vector* x) {
  double ans = 0.0;
  int offset = p*(p-1)/2;
  for(int n=0; n < offset; n++)
    ans += alphaij_logderiv(gsl_vector_get(x, n));
  for(int n=offset; n < x->size; n++) {
    double xn = gsl_vector_get(x, n);
    ans -= 0.5 * ( exp(xn) - xn ); /*1 deg. of freedom chisq. distrib. on log(x)*/
  }
  return ans;
}

double mcmclib_mcar_model_alphasigma_lpdf(void* in_p, gsl_vector* alphasigma) {
  mcmclib_mcar_model* p = (mcmclib_mcar_model*) in_p;

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
