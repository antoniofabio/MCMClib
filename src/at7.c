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
#include "at7.h"
#include "mvnorm.h"
#include "matrix.h"

/** \brief AT7 gamma values */
typedef struct {
  gsl_vector* beta; /**< current mixture weights*/
  gsl_vector** mu; /**< current mixture means*/
  gsl_matrix** Sigma; /**< current mixture variances*/

  mcmclib_mvnorm_lpdf** pik; /**< single mixture components densities*/
  mcmclib_mixnorm_lpdf* pi; /**< fitted mixture density*/

  gsl_matrix** qVariances; /**< local proposal variances*/
  mcmclib_mvnorm_lpdf** qdk; /**< proposal density components*/

  gsl_matrix* Sigma_eps; /**< positive-definiteness correction factor*/
  gsl_vector* scaling_factors; /**< region-specific scaling factors*/

  gsl_vector* tmpMean; /**< internal workspace memory*/
  gsl_vector* weights; /**< internal workspace memory*/
} at7_gamma;

/** \brief AT7 sufficient data */
typedef mcmclib_mixem_online* at7_suff;

#define AT7_GAMMA(p) ((at7_gamma*)((mcmclib_amh*) p)->mh->q->gamma)
#define AT7_SUFF(p) ((at7_suff*)((mcmclib_amh*) p)->suff)

static void at7_gamma_update_Sigma(at7_gamma* p) {
  const size_t K = p->beta->size;
  for(size_t k=0; k < K; k++) {
    mcmclib_matrix_addscale(p->qVariances[k],
			    p->Sigma_eps, p->Sigma[k],
			    gsl_vector_get(p->scaling_factors, k));
  }
}

/** alloc a new at7_gamma object. Input arguments are copied @internal */
static at7_gamma* at7_gamma_alloc(const gsl_vector* beta,
				  gsl_vector** mu,
				  gsl_matrix** Sigma) {
  at7_gamma* ans = (at7_gamma*) malloc(sizeof(at7_gamma));
  const size_t K = beta->size;
  const size_t dim = mu[0]->size;
  ans->beta = gsl_vector_alloc(K);
  gsl_vector_memcpy(ans->beta, beta);
  ans->mu = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  ans->Sigma = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  ans->pik = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  ans->tmpMean = gsl_vector_alloc(dim);
  ans->qVariances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  ans->qdk = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(size_t k=0; k < K; k++) {
    ans->mu[k] = gsl_vector_alloc(dim);
    gsl_vector_memcpy(ans->mu[k], mu[k]);
    ans->Sigma[k] = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(ans->Sigma[k], Sigma[k]);
    ans->pik[k] = mcmclib_mvnorm_lpdf_alloc(ans->mu[k], ans->Sigma[k]->data);
    ans->qVariances[k] = gsl_matrix_alloc(dim, dim);
    ans->qdk[k] = mcmclib_mvnorm_lpdf_alloc(ans->tmpMean, ans->qVariances[k]->data);
  }
  gsl_vector_memcpy(ans->beta, beta);
  ans->pi = mcmclib_mixnorm_lpdf_alloc(ans->beta, ans->pik);
  ans->Sigma_eps = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_identity(ans->Sigma_eps);
  gsl_matrix_scale(ans->Sigma_eps, 0.001);
  ans->scaling_factors = gsl_vector_alloc(K);
  gsl_vector_set_all(ans->scaling_factors, 2.38*2.38 / (double) dim);
  ans->weights = gsl_vector_alloc(K);
  gsl_vector_set_all(ans->weights, 0.0);
  at7_gamma_update_Sigma(ans);
  return ans;
}

/** frees an at7_gamma object @internal*/
static void at7_gamma_free(void* in_p) {
  if(!in_p) return;
  at7_gamma* p = (at7_gamma*) in_p;
  for(size_t k=0; k < p->beta->size; k++) {
    mcmclib_mvnorm_lpdf_free(p->pik[k]);
    gsl_vector_free(p->mu[k]);
    gsl_matrix_free(p->Sigma[k]);
    mcmclib_mvnorm_lpdf_free(p->qdk[k]);
    gsl_matrix_free(p->qVariances[k]);
  }
  mcmclib_mixnorm_lpdf_free(p->pi);
  gsl_matrix_free(p->Sigma_eps);
  gsl_vector_free(p->beta);
  gsl_vector_free(p->tmpMean);
  gsl_vector_free(p->weights);
  gsl_vector_free(p->scaling_factors);
  free(p->pik);
  free(p->qdk);
  free(p->mu);
  free(p->Sigma);
  free(p->qVariances);
  free(p);
}

static void at7_components_weights(at7_gamma* g, const gsl_vector* x) {
  mcmclib_mvnorm_lpdf** pik = g->pik;
  double sum=0.0;
  for(size_t k=0; k < g->weights->size; k++) {
    double nuk = gsl_vector_get(g->beta, k) * exp(mcmclib_mvnorm_lpdf_compute(pik[k], x));
    gsl_vector_set(g->weights, k, nuk);
    sum += nuk;
  }
  gsl_vector_scale(g->weights, 1.0 / sum);
}

/*sample a discrete value from the discrete distribution with probs. 'probs'*/
static size_t sample(gsl_rng* r, const gsl_vector* probs) {
  const size_t K = probs->size;
  double cum_sum = 0.0;
  const double who = gsl_rng_uniform(r);
  for(size_t which=0; which<K; which++) {
    if(who < (cum_sum += gsl_vector_get(probs, which)))
      return(which);
  }
  return(K-1);
}

static void at7_q_sample(mcmclib_mh_q* q, gsl_vector* x) {
  at7_gamma* gamma = (at7_gamma*) q->gamma;
  at7_components_weights(gamma, x);
  const size_t k = sample(q->r, gamma->weights);
  gsl_vector* oldx = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(oldx, x);
  mcmclib_mvnorm(q->r, gamma->qVariances[k], x);
  gsl_vector_add(x, oldx);
  gsl_vector_free(oldx);
}

static double at7_q_d(mcmclib_mh_q* q, gsl_vector* x, gsl_vector* y) {
  at7_gamma* gamma = (at7_gamma*) q->gamma;
  at7_components_weights(gamma, x);
  gsl_vector_memcpy(gamma->tmpMean, x);
  const size_t K = gamma->beta->size;
  double ans = 0.0;
  for(size_t k=0; k < K; k++) {
    double xy = exp(mcmclib_mvnorm_lpdf_compute(gamma->qdk[k], y));
    ans +=  xy * gsl_vector_get(gamma->weights, k);
  }
  return log(ans);
}

/** alloc a new AT7 sampler suff stats object
@param t0 burn-in length before starting adaptation
@returns a new AT7_suff object
*/
at7_suff* at7_suff_alloc(at7_gamma* g, size_t t0) {
  return (at7_suff*) mcmclib_mixem_online_alloc(g->mu, g->Sigma, g->beta, 0.5, t0);
}

/** free raptor_suff data*/
static void at7_suff_free(void* p) {
  mcmclib_mixem_online_free((mcmclib_mixem_online*) p);
}

/** Update suff stats of an AT7 chain*/
static void at7_suff_update(void* suff, const gsl_vector* x) {
  mcmclib_mixem_online* em = (mcmclib_mixem_online*) suff;
  mcmclib_mixem_online_update(em, x);
}

static void at7_gamma_update(void* suff, void* gamma) {
  suff = NULL; /*keep compiler quiet*/
  /*'suff' contains references to 'gamma' stuff*/
  at7_gamma_update_Sigma(gamma);
}

mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p logdistr, void* logdistr_data,
			       gsl_vector* x, size_t t0,
			       const gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat){
  at7_gamma* gamma = at7_gamma_alloc(beta_hat,
				     mu_hat,
				     Sigma_hat);
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, at7_q_sample, at7_q_d,
				       gamma, at7_gamma_free);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  at7_suff* suff = at7_suff_alloc(gamma, t0);
  mcmclib_amh* ans = mcmclib_amh_alloc(mh, t0, suff,
				       at7_suff_free,
				       at7_suff_update,
				       at7_gamma_update);
  return ans;
}

void mcmclib_at7_set_sf_all(mcmclib_amh* p, const double sf) {
  gsl_vector_set_all(AT7_GAMMA(p)->scaling_factors, sf);
}

void mcmclib_at7_set_sf(mcmclib_amh* p, const gsl_vector* sf) {
  gsl_vector_memcpy(AT7_GAMMA(p)->scaling_factors, sf);
}

/*here just for debugging purposes*/
#define PRINT_STATE(p)\
  fprintf(stderr, "w[0]: %f; means: mu[0]=%f, mu[1]=%f; variances: S[0]=%f, S[1]=%f\n",\
	  gsl_vector_get(AT7_GAMMA(p)->beta, 0),			\
	  gsl_vector_get(AT7_GAMMA(p)->mu[0], 0),			\
	  gsl_vector_get(AT7_GAMMA(p)->mu[1], 0),			\
	  gsl_matrix_get(AT7_GAMMA(p)->Sigma[0], 0, 0),			\
	  gsl_matrix_get(AT7_GAMMA(p)->Sigma[1], 0, 0))
