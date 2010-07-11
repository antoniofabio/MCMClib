/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <gsl/gsl_math.h>
#include "vector_stats.h"

void mcmclib_matrix_colmeans(const gsl_matrix* m, gsl_vector* out) {
  for(size_t i=0; i<m->size2; i++) {
    gsl_vector_const_view cv = gsl_matrix_const_column(m, i);
    gsl_vector_set(out, i, gsl_stats_mean(cv.vector.data, cv.vector.stride, cv.vector.size));
  }
}

void mcmclib_matrix_rowmeans(const gsl_matrix* m, gsl_vector* out) {
  for(size_t i=0; i<m->size1; i++) {
    gsl_vector_const_view rv = gsl_matrix_const_row(m, i);
    gsl_vector_set(out, i, gsl_stats_mean(rv.vector.data, rv.vector.stride, rv.vector.size));
  }
}

void mcmclib_matrix_covariance(const gsl_matrix* m, gsl_matrix* out) {
	const size_t d = m->size2;
	const size_t n = m->size1;
	gsl_matrix* mean = gsl_matrix_alloc(1, d);
	const gsl_matrix* row;

	gsl_vector_view mv = gsl_matrix_row(mean, 0);
	mcmclib_matrix_colmeans(m, &(mv.vector));

	gsl_matrix_set_zero(out);
	for(size_t i=0; i<n; i++) {
		gsl_matrix_const_view rv = gsl_matrix_const_submatrix (m, i, 0, 1, d);
		row = &(rv.matrix);
		gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, row, row, 1.0, out);
	}

	gsl_blas_dgemm (CblasTrans, CblasNoTrans, (double) -n, mean, mean, 1.0, out);

	gsl_matrix_scale(out, 1.0 / (double) n);

	gsl_matrix_free(mean);
}

void mcmclib_covariance_update(gsl_matrix* cov, gsl_vector* mean, size_t* n, const gsl_vector* x) {
  const size_t d = cov->size1;
  gsl_matrix_view colmean_view = gsl_matrix_view_array(mean->data, d, 1);
  gsl_matrix* colmean = &(colmean_view.matrix);
  gsl_matrix_const_view colx_view = gsl_matrix_const_view_array(x->data, d, 1);
  const gsl_matrix* colx = &(colx_view.matrix);

  if((*n) == 0) { /*this is the first call: do some cleanup*/
    gsl_vector_set_all(mean, 0.0);
    gsl_matrix_set_all(cov, 0.0);
  }

  /*update X %*% t(X) value:*/
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, (double) (*n), colmean, colmean, (double) (*n), cov);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, colx, colx, 1.0, cov);

  /*update mean value*/
  gsl_vector_scale(mean, (double) (*n));
  gsl_vector_add(mean, x);
  (*n)++;
  gsl_vector_scale(mean, 1.0 / (double) (*n));

  /*update covariance value*/
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, (double) -(*n), colmean, colmean, 1.0, cov);
  gsl_matrix_scale(cov, 1.0 / (double) (*n));
}

/*Pooled weighted variance*/
void mcmclib_pooled_variance(const double beta,
			     const gsl_vector** means,
			     const gsl_matrix** variances,
			     gsl_matrix* V) {
  const size_t dim = means[0]->size;
  gsl_matrix_memcpy(V, variances[0]);
  gsl_matrix_scale(V, beta);
  gsl_matrix* tmp = gsl_matrix_alloc(dim, dim);
  gsl_matrix_memcpy(tmp, variances[1]);
  gsl_matrix_scale(tmp, 1.0 - beta);
  gsl_matrix_add(V, tmp);

  gsl_matrix_view mu1v = gsl_matrix_view_array(means[0]->data, dim, 1);
  gsl_matrix* mu1 = &(mu1v.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mu1, mu1, 0.0, tmp);
  gsl_matrix_view mu2v = gsl_matrix_view_array(means[1]->data, dim, 1);
  gsl_matrix* mu2 = &(mu2v.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mu2, mu2, 1.0, tmp);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, mu1, mu2, 1.0, tmp);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, mu2, mu1, 1.0, tmp);
  gsl_matrix_scale(tmp, beta * (1-beta));
  gsl_matrix_add(V, tmp);

  gsl_matrix_free(tmp);
}

int mcmclib_vector_finite(gsl_vector* x) {
  for(size_t i=0; i<x->size; i++)
    if(!gsl_finite(gsl_vector_get(x, i)))
      return 0;
  return 1;
}
