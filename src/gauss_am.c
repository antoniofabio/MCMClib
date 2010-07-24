/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <gsl/gsl_blas.h>
#include <assert.h>
#include "gauss_am.h"
#include "gauss_mrw.h"

/**\brief Gaussian AM cumulated sufficient statistics and support data*/
typedef struct {
  size_t t; /**< number of iterations so far */
  gsl_vector* sum_x; /**< cumulated sum of xs*/
  gsl_matrix* sum_xx; /**< cumulated sum of xxs*/
  gsl_matrix* Sigma_eps; /**< pos. definiteness cov. correction additive constant*/
  gsl_matrix* Sigma_zero; /**< starting proposal covariance matrix*/
  double sf; /**< scaling factor*/
} gauss_am_suff;

/** free extra AM data*/
static void gauss_am_suff_free(void* p) {
  if(!p) return;
  gauss_am_suff* s = (gauss_am_suff*) p;
  gsl_matrix_free(s->Sigma_eps);
  gsl_matrix_free(s->Sigma_zero);
  gsl_vector_free(s->sum_x);
  gsl_matrix_free(s->sum_xx);
  free(s);
}

static void gauss_am_update_suff(void* in_s, const gsl_vector* x) {
  gauss_am_suff* s = (gauss_am_suff*) in_s;
  s->t++;
  gsl_vector_add(s->sum_x, x);
  gsl_matrix_const_view x_cv = gsl_matrix_const_view_array(x->data, x->size, 1);
  const gsl_matrix* x_cm = &(x_cv.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, x_cm, x_cm, 1.0, s->sum_xx);
}

/** AM gamma update function \internal */
static void gauss_am_update_gamma(void* in_suff, void* in_gamma) {
  gauss_am_suff* s = (gauss_am_suff*) in_suff;
  gsl_matrix* Sigma = (gsl_matrix*) in_gamma;
  const size_t t = s->t;

  const size_t d = Sigma->size1;
  gsl_vector* mean = gsl_vector_alloc(d);
  gsl_vector_memcpy(mean, s->sum_x);
  gsl_vector_scale(mean, 1.0 / (double) t);
  gsl_matrix_view mean_cv = gsl_matrix_view_array(mean->data, d, 1);
  gsl_matrix* mean_cm = &(mean_cv.matrix);

  gsl_matrix_memcpy(Sigma, s->sum_xx);
  gsl_matrix_scale(Sigma, 1.0 / (double) t);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, mean_cm, mean_cm, 1.0, Sigma);
  gsl_matrix_add(Sigma, s->Sigma_eps);
  gsl_matrix_scale(Sigma, s->sf);
  gsl_vector_free(mean);
}

mcmclib_amh* mcmclib_gauss_am_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* logdistr_data,
				    gsl_vector* start_x,
				    const gsl_matrix* sigma_zero, size_t t0) {
  gauss_am_suff* suff = (gauss_am_suff*) malloc(sizeof(gauss_am_suff));
  assert(sigma_zero->size1 == sigma_zero->size2);
  assert(start_x->size == sigma_zero->size1);
  const size_t d = start_x->size;
  suff->Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_memcpy(suff->Sigma_zero, sigma_zero);
  suff->t = 0;
  suff->sum_x = gsl_vector_alloc(d);
  gsl_vector_set_zero(suff->sum_x);
  suff->sum_xx = gsl_matrix_alloc(d, d);
  gsl_matrix_set_zero(suff->sum_xx);
  suff->sf = (2.38 * 2.38) / (double) d;

  suff->Sigma_eps = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(suff->Sigma_eps);
  gsl_matrix_scale(suff->Sigma_eps, 0.001);

  return mcmclib_amh_alloc(mcmclib_gauss_mrw_alloc(r, logdistr, logdistr_data,
						   start_x, sigma_zero),
			   t0, suff,
			   gauss_am_suff_free, gauss_am_update_suff, gauss_am_update_gamma);
}

void mcmclib_gauss_am_set_sf(mcmclib_amh* p, double sf) {
  assert(sf > 0.0);
  ((gauss_am_suff*) (p->suff))->sf = sf;
}

void mcmclib_gauss_am_set_eps(mcmclib_amh* p, double eps) {
  assert(eps > 0.0);
  gauss_am_suff* s = (gauss_am_suff*) p->suff;
  gsl_matrix_set_identity(s->Sigma_eps);
  gsl_matrix_scale(s->Sigma_eps, eps);
}
