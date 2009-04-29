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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "spatial.h"

mcmclib_spatial_lpdf* mcmclib_spatial_lpdf_alloc(gsl_vector* mu,
						 gsl_vector* rho,
						 gsl_vector* sigma,
						 gsl_vector* tausq,
						 gsl_matrix* D) {
  mcmclib_spatial_lpdf* a = (mcmclib_spatial_lpdf*) malloc(sizeof(mcmclib_spatial_lpdf));
  a->mu = mu;
  assert(mu->size == D->size1);
  a->rho = rho;
  assert(rho->size == 1);
  a->sigma = sigma;
  assert(sigma->size == 1);
  a->tausq = tausq;
  assert(tausq->size == 1);
  a->D = D;
  assert(D->size1 == D->size2);
  a->Sigma = gsl_matrix_alloc(mu->size, mu->size);
  a->norm = mcmclib_mvnorm_lpdf_alloc(mu, a->Sigma->data);
  return a;
}

void mcmclib_spatial_lpdf_free(mcmclib_spatial_lpdf* p) {
  gsl_matrix_free(p->Sigma);
  mcmclib_mvnorm_lpdf_free(p->norm);
  free(p);
}

double mcmclib_spatial_cov_exponential(double d, double rho, double sigma, double tausq) {
  double a = (sigma - tausq) * (1.0 - exp(- d / rho));
  if (d>0)
    a += tausq;
  return a;
}

void mcmclib_spatial_distances(gsl_matrix* D, gsl_matrix* xy) {
  int n = xy->size1;
  int d = xy->size2;
  gsl_matrix_set_all(D, 0.0);
  gsl_vector* x = gsl_vector_alloc(d);
  for(int i=0; i<n; i++) for(int j=i+1; j<n; j++) {
      gsl_vector_view ri_v = gsl_matrix_row(xy, i);
      gsl_vector* ri = &(ri_v.vector);
      gsl_vector_memcpy(x, ri);
      gsl_vector_view rj_v = gsl_matrix_row(xy, j);
      gsl_vector* rj = &(rj_v.vector);
      gsl_vector_sub(x, rj);
      gsl_vector_mul(x, x);
      double dist = 0.0;
      for(int k = 0; k<d; k++)
	dist += gsl_vector_get(x, k);
      dist = sqrt(dist);
      gsl_matrix_set(D, i, j, dist);
      gsl_matrix_set(D, j, i, dist);
  }
  gsl_vector_free(x);
}

void mcmclib_spatial_set_xy(mcmclib_spatial_lpdf* p, gsl_matrix* xy) {
  mcmclib_spatial_distances(p->D, xy);
}

/*check if the real, simm. matrix A is pos. def.*/
int matrix_posDef(gsl_matrix* A) {
  int n = A->size1;
  gsl_eigen_symm_workspace* work = gsl_eigen_symm_alloc(n);
  gsl_matrix* A1 = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(A1, A);
  gsl_vector* v = gsl_vector_alloc(n);
  gsl_eigen_symm(A1, v, work);
  gsl_eigen_symm_free(work);
  gsl_matrix_free(A1);
  int ans = gsl_vector_ispos(v);
  gsl_vector_free(v);
  return ans;
}

double mcmclib_spatial_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_spatial_lpdf* p = (mcmclib_spatial_lpdf*) in_p;
  int n = x->size;
  gsl_matrix* D = p->D;
  assert(n == D->size1);
  double rho = gsl_vector_get(p->rho, 0);
  double sigma = gsl_vector_get(p->sigma, 0);
  double tausq = gsl_vector_get(p->tausq, 0);
  for(int i=0; i<n; i++)
    gsl_matrix_set(p->Sigma, i, i, sigma);
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++) {
      double gammaij = mcmclib_spatial_cov_exponential(gsl_matrix_get(D, i, j), rho, sigma, tausq);
      gsl_matrix_set(p->Sigma, i, j, sigma - gammaij);
      gsl_matrix_set(p->Sigma, j, i, sigma - gammaij);
    }
  for(int i=0; i<n; i++)
    if(gsl_matrix_get(p->Sigma, i, i) < 1e-6)
      return log(0.0);
  if(!matrix_posDef(p->Sigma))
    return log(0.0);
  return mcmclib_mvnorm_lpdf_compute(p->norm, x);
}
