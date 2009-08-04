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
  double within = 0.0;
  for(int k=0; k<K; k++) {
    gsl_matrix_memcpy(work, g->Sigma_hat[k]);
    if(mcmclib_cholesky_decomp(work) == GSL_SUCCESS)
      within += exp(mcmclib_matrix_logtrace(work));
  }
  gsl_vector* mean = gsl_vector_alloc(d);
  int n=0;
  for(int k=0; k<K; k++)
    mcmclib_covariance_update(work, mean, &n, g->mu_hat[k]);
  double between = 0.0;
  if(mcmclib_cholesky_decomp(work) == GSL_SUCCESS)
    between = exp(mcmclib_matrix_logtrace(work));
  gsl_matrix_free(work);
  gsl_vector_free(mean);
  if(between == within)
    return 0.5;
  return between / (between + within);
}

double mcmclib_raptor_alpha_identity_fun(void* ignore, mcmclib_raptor_gamma* g) {
  return mcmclib_raptor_alpha_star_fun(g);
}

void mcmclib_raptor_set_alpha_fun_identity(mcmclib_amh* p) {
  mcmclib_raptor_set_alpha_fun(p, NULL, mcmclib_raptor_alpha_identity_fun);
}
