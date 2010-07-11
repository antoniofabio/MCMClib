/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
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
#include "matrix.h"
#include "lpdf_iwishart.h"

mcmclib_iwishart_lpdf* mcmclib_iwishart_lpdf_alloc(const gsl_matrix* Psi, const unsigned int m) {
  mcmclib_iwishart_lpdf* a = (mcmclib_iwishart_lpdf*) malloc(sizeof(mcmclib_iwishart_lpdf));
  assert(Psi->size1 == Psi->size2);
  const size_t n = Psi->size1;
  a->Psi = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(a->Psi, Psi);
  gsl_matrix* PsiChol = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(PsiChol, Psi);
  gsl_linalg_cholesky_decomp(PsiChol);
  a->PsiDet = 0.0;
  for(size_t i=0; i<n; i++)
    a->PsiDet += log(gsl_matrix_get(PsiChol, i, i));
  gsl_matrix_free(PsiChol);
  a->PsiDet *= (double) m;
  a->m = m;

  a->PsiX = gsl_matrix_alloc(n, n);
  gsl_matrix_set_zero(a->PsiX);
  a->X1 = gsl_matrix_alloc(n, n);
  return a;
}

void mcmclib_iwishart_lpdf_free(mcmclib_iwishart_lpdf* p) {
  if(!p) return;
  gsl_matrix_free(p->X1);
  gsl_matrix_free(p->PsiX);
  gsl_matrix_free(p->Psi);
  free(p);
}

static double trace(const gsl_matrix* A) {
  double ans = 0.0;
  for(size_t i=0; i < A->size1; i++)
    ans += gsl_matrix_get(A, i, i);
  return ans;
}

double mcmclib_iwishart_lpdf_compute(void* in_p, const gsl_vector* x) {
  mcmclib_iwishart_lpdf* p = (mcmclib_iwishart_lpdf*) in_p;
  const size_t n = p->Psi->size1;
  gsl_matrix_const_view X_v = gsl_matrix_const_view_vector(x, n, n);
  const gsl_matrix* X = &(X_v.matrix);
  gsl_matrix* X1 = p->X1;
  gsl_matrix_memcpy(X1, X);
  if(mcmclib_cholesky_decomp(X1) != GSL_SUCCESS)
    return log(0.0);
  double Xdet2 = 0.0;
  for(size_t i=0; i<n; i++)
    Xdet2 += log(gsl_matrix_get(X1, i, i));
  gsl_linalg_cholesky_invert(X1);
  gsl_matrix* PsiX = p->PsiX;
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, p->Psi, X1, 0.0, PsiX);
  return p->PsiDet - ((double) p->m + (double) n + 1.0) * Xdet2 -0.5 * trace(PsiX);
}
