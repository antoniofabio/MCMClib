/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include "raptor.h"
#include "matrix.h"
#include "vector_stats.h"

double mcmclib_raptor_alpha_default_fun(void* data, mcmclib_raptor_gamma* p) {
  mcmclib_rapt_gamma* g = (mcmclib_rapt_gamma*) data;
  return gsl_matrix_get(g->lambda, 0, g->lambda->size2 - 1);
}

double mcmclib_raptor_alpha_star_fun(mcmclib_raptor_gamma* g) {
  int K = g->beta_hat->size;
  int d = g->mu_hat[0]->size;
  gsl_matrix* work = gsl_matrix_alloc(d, d);
  gsl_matrix* W = gsl_matrix_alloc(d, d);
  gsl_matrix_set_zero(W);
  for(int k=0; k<K; k++) {
    gsl_matrix_memcpy(work, g->Sigma_hat[k]);
    gsl_matrix_scale(work, gsl_vector_get(g->beta_hat, k));
    gsl_matrix_add(W, work);
  }
  double lwithin = log(0.0);
  gsl_matrix_free(work);
  if(mcmclib_cholesky_decomp(W) == GSL_SUCCESS) {
    lwithin = mcmclib_matrix_logtrace(W);
  } else {
    gsl_matrix_free(W);
    return 1.0;
  }
  gsl_matrix_free(W);
  gsl_matrix* B = gsl_matrix_alloc(d, d);
  gsl_matrix_set_zero(B);
  for(int k=0; k<K; k++) {
    gsl_matrix_view Mkv = gsl_matrix_view_vector(g->mu_hat[k], d, 1);
    gsl_matrix* Mk = &(Mkv.matrix);
    double wk = gsl_vector_get(g->beta_hat, k);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, wk * (1.0 - wk),
		    Mk, Mk, 1.0, B);
  }
  double lbetween = 0.0;
  if(mcmclib_cholesky_decomp(B) == GSL_SUCCESS)
    lbetween = mcmclib_matrix_logtrace(B);
  gsl_matrix_free(B);
  if(lbetween == lwithin)
    return 0.5;
  return exp(lbetween) / (exp(lbetween) + exp(lwithin));
}

double mcmclib_raptor_alpha_identity_fun(void* ignore, mcmclib_raptor_gamma* g) {
  return mcmclib_raptor_alpha_star_fun(g);
}

void mcmclib_raptor_set_alpha_fun_identity(mcmclib_amh* p) {
  mcmclib_raptor_set_alpha_fun(p, NULL, mcmclib_raptor_alpha_identity_fun);
}
