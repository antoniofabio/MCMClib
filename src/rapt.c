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
#include "rapt.h"
#include "mvnorm.h"
#include "vector_stats.h"

/**\brief RAPT sufficient statistics*/
typedef struct {
  size_t t; /**current iteration number*/
  gsl_vector** means; /**< array of regions means*/
  gsl_matrix** variances; /**< array of regions variances*/
  gsl_vector* global_mean;
  gsl_matrix* global_variance;
  gsl_vector* n; /**< number of visits in each region*/

  /*internal data*/
  gsl_matrix* Sigma_eps; /**< additive perturbation factor for variances updating*/
  gsl_vector* workspace; /**< utility workspace memory*/
  double scaling_factor_local; /**< local proposal variance scaling factor*/
  double scaling_factor_global; /**< global proposal variance scaling factor*/
  
  region_fun_t which_region_fun;
  void* which_region_data;
} rapt_suff;

static void rapt_suff_set_correction_factor(rapt_suff* p, const double eps) {
  gsl_matrix_set_identity(p->Sigma_eps);
  gsl_matrix_scale(p->Sigma_eps, eps);
}

static void rapt_suff_set_scaling_factors(rapt_suff* p, const double local, const double global) {
  p->scaling_factor_local = local;
  p->scaling_factor_global = global;
}

static rapt_suff* rapt_suff_alloc(mcmclib_rapt_gamma* gamma) {
  rapt_suff* a = (rapt_suff*) malloc(sizeof(rapt_suff));
  const size_t K = gamma->K;
  const size_t dim = gamma->sigma_whole->size1;
  a->t = 0;
  a->means = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  a->variances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(size_t k=0; k<K; k++) {
    a->means[k] = gsl_vector_alloc(dim);
    a->variances[k] = gsl_matrix_alloc(dim, dim);
  }
  a->global_mean = gsl_vector_alloc(dim);
  gsl_vector_set_all(a->global_mean, 0.0);
  a->global_variance = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_all(a->global_variance, 0.0);
  a->n = gsl_vector_alloc(K);
  gsl_vector_set_all(a->n, 0.0);
  a->Sigma_eps = gsl_matrix_alloc(dim, dim);

  rapt_suff_set_correction_factor(a, 0.001);
  double sf = 2.38 * 2.38 / ((double) dim);
  rapt_suff_set_scaling_factors(a, sf, sf);

  a->which_region_fun = gamma->which_region;
  a->which_region_data = gamma->which_region_data;
  return a;
}

static void rapt_suff_free(void* in_p) {
  rapt_suff* p = (rapt_suff*) in_p;
  size_t K = p->n->size;
  gsl_vector_free(p->n);
  for(size_t k=0; k< K; k++) {
    gsl_vector_free(p->means[k]);
    gsl_matrix_free(p->variances[k]);
  }
  gsl_vector_free(p->global_mean);
  gsl_matrix_free(p->global_variance);
  gsl_matrix_free(p->Sigma_eps);
  free(p->means);
  free(p->variances);

  free(p);
}

static void rapt_suff_update(void* in_suff, const gsl_vector* x) {
  rapt_suff* s = (rapt_suff*) in_suff;
  s->t ++;
  const size_t k = s->which_region_fun(s->which_region_data, x);
  size_t fake_n = (size_t) gsl_vector_get(s->n, k);
  mcmclib_covariance_update(s->variances[k], s->means[k], &fake_n, x);
  fake_n = (size_t) (s->t - 1);
  mcmclib_covariance_update(s->global_variance, s->global_mean, &fake_n, x);
  gsl_vector_set(s->n, k, gsl_vector_get(s->n, k) + 1);
}

static void rapt_gamma_update(void* suff, void* gamma) {
  rapt_suff* s = (rapt_suff*) suff;
  mcmclib_rapt_gamma* g = (mcmclib_rapt_gamma*) gamma;
  mcmclib_rapt_q_update_proposals_custom(g,
					 s->variances, s->global_variance,
					 s->Sigma_eps,
					 s->scaling_factor_local,
					 s->scaling_factor_global);
}

mcmclib_amh* mcmclib_rapt_alloc(gsl_rng* r,
				distrfun_p logdistr, void* logdistr_data,
				gsl_vector* x,
				const size_t t0,
				const gsl_matrix* sigma_whole,
				size_t K,
				gsl_matrix** sigma_local,
				region_fun_t which_region,
				void* which_region_data) {
  mcmclib_mh_q* q = mcmclib_rapt_q_alloc(r,
					 sigma_whole, K, sigma_local,
					 which_region, which_region_data);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  rapt_suff* suff = rapt_suff_alloc(q->gamma);

  return mcmclib_amh_alloc(mh, t0, suff, rapt_suff_free, rapt_suff_update,
			   rapt_gamma_update);
}

void mcmclib_rapt_set_correction_factor(mcmclib_amh* p, const double eps) {
  rapt_suff* s = (rapt_suff*) p->suff;
  rapt_suff_set_correction_factor(s, eps);
}

void mcmclib_rapt_set_scaling_factors(mcmclib_amh* p, const double local, const double global) {
  rapt_suff* s = (rapt_suff*) p->suff;
  rapt_suff_set_scaling_factors(s, local, global);
}
