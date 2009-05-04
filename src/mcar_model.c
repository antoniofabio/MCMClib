/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_vector.h>
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

static int is_sorted(gsl_vector* v) {
  double m = gsl_vector_get(v, 0);
  for(int i=1; i<v->size; i++) {
    double n = gsl_vector_get(v, i);
    if(n > m)
      return 0;
    m = n;
  }
  return 1;
}

double mcmclib_mcar_model_alpha12sigma_lpdf(void* in_p, gsl_vector* alpha12sigma) {
  mcmclib_mcar_model* p = (mcmclib_mcar_model*) in_p;
  int n = p->lpdf->p;
  gsl_vector_view s_v = gsl_vector_subvector(alpha12sigma, n*(n-1), n);
  if(!is_sorted(&s_v.vector))
    return log(0.0);
  gsl_vector* tmp = gsl_vector_alloc(alpha12sigma->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha12sigma);
  gsl_vector_memcpy(p->lpdf->alpha12sigma, alpha12sigma);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha12sigma, tmp);
  gsl_vector_free(tmp);
  return ans;
}

static void mprint(gsl_matrix* A) {
  int n = A->size1;
  int p = A->size2;
  for(int i=0; i<n; i++) {
    for(int j=0; j<p; j++)
      printf("%.3f, ", gsl_matrix_get(A, i, j));
    printf("\n");
  }
}

static double iwishart_alphasigma(mcmclib_iwishart_lpdf* p, gsl_vector* as) {
  int n = p->Psi->size1;
  gsl_vector_view s_v = gsl_vector_subvector(as, n*(n-1)/2, n);
  if(!is_sorted(&s_v.vector))
    return log(0.0);
  gsl_matrix* Gamma = gsl_matrix_alloc(n, n);
  mcmclib_Givens_representation(Gamma, as);
  gsl_vector_view gamma_v = gsl_vector_view_array(Gamma->data, n*n);
  double ans = mcmclib_iwishart_lpdf_compute(p, &gamma_v.vector);
  gsl_matrix_free(Gamma);
  return ans;
}

double mcmclib_mcar_model_alphasigma_lpdf(void* in_p, gsl_vector* alphasigma) {
  mcmclib_mcar_model* p = (mcmclib_mcar_model*) in_p;

  double ans = iwishart_alphasigma(p->w, alphasigma);
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
