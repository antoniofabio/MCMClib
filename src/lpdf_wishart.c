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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "lpdf_wishart.h"

mcmclib_wishart_lpdf* mcmclib_wishart_lpdf_alloc(gsl_matrix* V, int p) {
  mcmclib_wishart_lpdf* a = (mcmclib_wishart_lpdf*) malloc(sizeof(mcmclib_wishart_lpdf));
  assert(V->size1 == V->size2);
  int n = V->size1;
  a->V = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(a->V, V);
  gsl_linalg_cholesky_decomp(a->V);
  gsl_linalg_cholesky_invert(a->V);
  a->p = p;

  a->Vx = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(a->Vx);
  a->X1 = gsl_matrix_alloc(n, n);
  return a;
}

void mcmclib_wishart_lpdf_free(mcmclib_wishart_lpdf* p) {
  gsl_matrix_free(p->X1);
  gsl_matrix_free(p->Vx);
  gsl_matrix_free(p->V);
  free(p);
}

double mcmclib_wishart_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_wishart_lpdf* p = (mcmclib_wishart_lpdf*) in_p;
  int n = p->V->size1;
  gsl_matrix_view X_v = gsl_matrix_view_vector(x, n, n);
  gsl_matrix* X = &(X_v.matrix);
  gsl_matrix* Vx = p->Vx;
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, p->V, X, 0.0, Vx);
  double ans = 0.0;
  for(int i=0; i<n; i++)
    ans += gsl_matrix_get(Vx, i, i);
  ans *= -0.5;
  gsl_matrix* X1 = p->X1;
  gsl_matrix_memcpy(X1, X);
  gsl_linalg_cholesky_decomp(X1);
  double ldet_X = 0.0;
  for(int i=0; i<n; i++)
    ldet_X += log(gsl_matrix_get(X1, i, i));
  ans += (n - p->p - 1.0) * ldet_X;
  return ans;
}
