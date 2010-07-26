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
#include "raptor.h"
#include "region_mixnorm.h"
#include "mvnorm.h"

/** \brief RAPTOR sampler gamma values */
typedef struct {
  gsl_vector* beta_hat; /**< current mixture weights estimates*/
  gsl_vector** mu_hat; /**< current mixture means estimates*/
  gsl_matrix** Sigma_hat; /**< current mixture variances estimates*/

  mcmclib_mvnorm_lpdf** pik_hat; /**< single mixture components densities*/
  mcmclib_mixnorm_lpdf* pi_hat; /**< mixture density*/
} raptor_gamma;

/** alloc a new raptor_gamma object. Input arguments are copied @internal */
static raptor_gamma* raptor_gamma_alloc(const gsl_vector* beta_hat,
					gsl_vector** mu_hat,
					gsl_matrix** Sigma_hat) {
  raptor_gamma* a = (raptor_gamma*) malloc(sizeof(raptor_gamma));
  a->beta_hat = gsl_vector_alloc(beta_hat->size);
  gsl_vector_memcpy(a->beta_hat, beta_hat);
  size_t K = beta_hat->size;
  size_t d = mu_hat[0]->size;
  a->mu_hat = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  a->Sigma_hat = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(size_t k=0; k < K; k++) {
    a->mu_hat[k] = gsl_vector_alloc(d);
    gsl_vector_memcpy(a->mu_hat[k], mu_hat[k]);
    a->Sigma_hat[k] = gsl_matrix_alloc(d, d);
    gsl_matrix_memcpy(a->Sigma_hat[k], Sigma_hat[k]);
  }

  a->pik_hat = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(size_t k=0; k<K; k++)
    a->pik_hat[k] = mcmclib_mvnorm_lpdf_alloc(a->mu_hat[k], a->Sigma_hat[k]->data);

  a->pi_hat = mcmclib_mixnorm_lpdf_alloc(a->beta_hat, a->pik_hat);
  return a;
}

/** frees a raptor_gamma object @internal*/
static void raptor_gamma_free(raptor_gamma* p) {
  mcmclib_mixnorm_lpdf_free(p->pi_hat);
  for(size_t k=0; k < p->beta_hat->size; k++) {
    mcmclib_mvnorm_lpdf_free(p->pik_hat[k]);
    gsl_vector_free(p->mu_hat[k]);
    gsl_matrix_free(p->Sigma_hat[k]);
  }
  free(p->mu_hat);
  free(p->Sigma_hat);
  gsl_vector_free(p->beta_hat);
  free(p->pik_hat);
  free(p);
}

/** \brief RAPTOR sufficient data */
typedef struct {
  mcmclib_mixem_online* em; /**< online-EM mixture fitter*/

  gsl_matrix* Sigma_eps;
  double scaling_factor_local;
  double scaling_factor_global;
} raptor_suff;

/** alloc a new RAPTOR sampler suff stats object
@param t0 burn-in length before starting adaptation
@returns a new raptor_suff object
*/
static raptor_suff* raptor_suff_alloc(raptor_gamma* g, size_t t0) {
  raptor_suff* a = (raptor_suff*) malloc(sizeof(raptor_suff));
  a->em = mcmclib_mixem_online_alloc(g->mu_hat, g->Sigma_hat, g->beta_hat, 0.6, t0);
  size_t d = g->mu_hat[0]->size;
  a->Sigma_eps = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(a->Sigma_eps);
  gsl_matrix_scale(a->Sigma_eps, 0.001);
  return a;
}

/** free raptor_suff data*/
static void raptor_suff_free(void* in_p) {
  raptor_suff* p = (raptor_suff*) in_p;
  mcmclib_mixem_online_free(p->em);
  gsl_matrix_free(p->Sigma_eps);
  free(p);
}

/** Update suff stats of a RAPTOR chain*/
static void raptor_suff_update(void* p, const gsl_vector* x) {
  raptor_suff* s = (raptor_suff*) p;
  mcmclib_mixem_online* em = s->em;
  mcmclib_mixem_online_update(em, x);  
}

#define RAPT_GAMMA(p) ((mcmclib_rapt_gamma*) (p))
#define RAPTOR_GAMMA(p) ((raptor_gamma*) RAPT_GAMMA(p)->which_region_data)

static void raptor_gamma_update(void* in_s, void* in_g) {
  raptor_suff* s = (raptor_suff*) in_s;
  mcmclib_rapt_q_update_proposals_custom(RAPT_GAMMA(in_g),
					 s->em->Sigma, s->em->Sigma_global,
					 s->Sigma_eps,
					 s->scaling_factor_local,
					 s->scaling_factor_global);
}

static size_t which_region_fun(void* in_g, const gsl_vector* x) {
  return mcmclib_region_mixnorm_compute(x, ((raptor_gamma*) in_g)->pi_hat);
}

mcmclib_amh* mcmclib_raptor_alloc(gsl_rng* r,
				  distrfun_p logdistr, void* logdistr_data,
				  gsl_vector* x, const size_t t0, gsl_matrix* Sigma_zero,
				  const gsl_vector* beta_hat,
				  gsl_vector** mu_hat,
				  gsl_matrix** Sigma_hat){
  size_t K = beta_hat->size;

  raptor_gamma* gamma = raptor_gamma_alloc(beta_hat,
					   mu_hat,
					   Sigma_hat);
  mcmclib_mh_q* q = mcmclib_rapt_q_alloc(r, Sigma_zero, K, Sigma_hat,
					 which_region_fun, gamma,
					 (free_fun_t) raptor_gamma_free);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  raptor_suff* suff = raptor_suff_alloc(gamma, t0);
  mcmclib_amh* ans = mcmclib_amh_alloc(mh, t0, suff, raptor_suff_free,
				       raptor_suff_update, raptor_gamma_update);

  mcmclib_raptor_set_sf(ans, 2.38*2.38 / (double) x->size);
  return ans;
}

void mcmclib_raptor_set_sf(mcmclib_amh* p, double sf) {
  mcmclib_raptor_set_sf_global(p, sf);
  mcmclib_raptor_set_sf_local(p, sf);
}

void mcmclib_raptor_set_sf_global(mcmclib_amh* p, double sf) {
  raptor_suff* s = (raptor_suff*) p->suff;
  s->scaling_factor_global = sf;
}

void mcmclib_raptor_set_sf_local(mcmclib_amh* p, double sf) {
  raptor_suff* s = (raptor_suff*) p->suff;
  s->scaling_factor_local = sf;
}

void mcmclib_raptor_set_alpha(mcmclib_amh* p, double alpha) {
  mcmclib_rapt_gamma_set_alpha((mcmclib_rapt_gamma*) p->mh->q->gamma, alpha);
}
