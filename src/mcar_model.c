/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "mcar_model.h"

mcmclib_mcar_model* mcmclib_mcar_model_alloc(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e) {
  mcmclib_mcar_model* a = (mcmclib_mcar_model*) malloc(sizeof(mcmclib_mcar_model));
  a->lpdf = m;
  a->e = e;
  return a;
}

void mcmclib_mcar_model_free(mcmclib_mcar_model* p) {
  free(p);
}

double mcmclib_mcar_model_alpha1_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha1) {
  for(int i=0; i < alpha1->size; i++) {
    double alpha1i = gsl_vector_get(alpha1, i);
    if((alpha1i <= -M_PI/2.0) || (alpha1i >= M_PI/2.0))
      return log(0.0);
  }
  gsl_vector* tmp = gsl_vector_alloc(alpha1->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha1);
  gsl_vector_memcpy(p->lpdf->alpha1, alpha1);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha1, tmp);
  gsl_vector_free(tmp);
  return ans;
}

double mcmclib_mcar_model_alpha2_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha2) {
  for(int i=0; i < alpha2->size; i++) {
    double alpha2i = gsl_vector_get(alpha2, i);
    if((alpha2i <= -M_PI/2.0) || (alpha2i >= M_PI/2.0))
      return log(0.0);
  }
  gsl_vector* tmp = gsl_vector_alloc(alpha2->size);
  gsl_vector_memcpy(tmp, p->lpdf->alpha2);
  gsl_vector_memcpy(p->lpdf->alpha2, alpha2);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->alpha2, tmp);
  gsl_vector_free(tmp);
  return ans;
}

double mcmclib_mcar_model_sigma_lpdf(mcmclib_mcar_model* p, gsl_vector* sigma) {
  double sigma0 = gsl_vector_get(sigma, 0);
  if((sigma0 <= 0.0) || (sigma0 >= 1.0))
    return log(0.0);
  gsl_vector* tmp = gsl_vector_alloc(sigma->size);
  gsl_vector_memcpy(tmp, p->lpdf->sigma);
  gsl_vector_memcpy(p->lpdf->sigma, sigma);
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_vector_memcpy(p->lpdf->sigma, tmp);
  gsl_vector_free(tmp);
  return ans;
}

double mcmclib_mcar_model_Gamma_lpdf(mcmclib_mcar_model* p, gsl_vector* gamma) {
  int P = p->lpdf->p;
  gsl_matrix* tmp = gsl_matrix_alloc(P, P);
  gsl_matrix_view gamma_v = gsl_matrix_view_vector(gamma, P, P);
  gsl_matrix_memcpy(tmp, p->lpdf->Gamma);
  gsl_matrix_memcpy(p->lpdf->Gamma, &(gamma_v.matrix));
  double ans = mcmclib_mcar_tilde_lpdf_compute(p->lpdf, p->e);
  gsl_matrix_memcpy(p->lpdf->Gamma, tmp);
  gsl_matrix_free(tmp);
  return ans;
}
