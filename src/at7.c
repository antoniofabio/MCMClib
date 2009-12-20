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
#include "at7.h"
#include "mvnorm.h"
#include "matrix.h"

#define AT7_GAMMA(p) ((mcmclib_at7_gamma*)((mcmclib_amh*) p)->mh->q->gamma)
#define AT7_SUFF(p) ((mcmclib_at7_suff*)((mcmclib_amh*) p)->suff)

void at7_gamma_update_Sigma(mcmclib_at7_gamma* p) {
  int K = p->beta->size;
  for(int k=0; k < K; k++) {
    mcmclib_matrix_addscale(p->qVariances[k],
			    p->Sigma_eps, p->Sigma[k],
			    gsl_vector_get(p->scaling_factors, k));
  }
}

mcmclib_at7_gamma* mcmclib_at7_gamma_alloc(const gsl_vector* beta,
					   gsl_vector** mu,
					   gsl_matrix** Sigma) {
  mcmclib_at7_gamma* ans = (mcmclib_at7_gamma*) malloc(sizeof(mcmclib_at7_gamma));
  int K = beta->size;
  int dim = mu[0]->size;
  ans->beta = gsl_vector_alloc(K);
  gsl_vector_memcpy(ans->beta, beta);
  ans->mu = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  ans->Sigma = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  ans->pik = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  ans->tmpMean = gsl_vector_alloc(dim);
  ans->qVariances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  ans->qdk = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k < beta->size; k++) {
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

void mcmclib_at7_gamma_free(void* in_p) {
  mcmclib_at7_gamma* p = (mcmclib_at7_gamma*) in_p;
  for(int k=0; k < p->beta->size; k++) {
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

static void at7_components_weights(mcmclib_at7_gamma* g, gsl_vector* x) {
  mcmclib_mvnorm_lpdf** pik = g->pik;
  double sum=0.0;
  for(int k=0; k < g->weights->size; k++) {
    double nuk = gsl_vector_get(g->beta, k) * exp(mcmclib_mvnorm_lpdf_compute(pik[k], x));
    gsl_vector_set(g->weights, k, nuk);
    sum += nuk;
  }
  gsl_vector_scale(g->weights, 1.0 / sum);
}

/*sample a discrete value from the discrete distribution with probs. 'probs'*/
static int sample(gsl_rng* r, gsl_vector* probs) {
  int K = probs->size;
  double cum_sum = 0.0;
  double who = gsl_rng_uniform(r);
  for(int which=0; which<K; which++) {
    if(who < (cum_sum += gsl_vector_get(probs, which)))
      return(which);
  }
  return(K-1);
}

static void at7_q_sample(mcmclib_mh_q* q, gsl_vector* x) {
  mcmclib_at7_gamma* gamma = (mcmclib_at7_gamma*) q->gamma;
  at7_components_weights(gamma, x);
  int k = sample(q->r, gamma->weights);
  gsl_vector* oldx = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(oldx, x);
  mcmclib_mvnorm(q->r, gamma->qVariances[k], x);
  gsl_vector_add(x, oldx);
  gsl_vector_free(oldx);
}

static double at7_q_d(mcmclib_mh_q* q, gsl_vector* x, gsl_vector* y) {
  mcmclib_at7_gamma* gamma = (mcmclib_at7_gamma*) q->gamma;
  at7_components_weights(gamma, x);
  gsl_vector_memcpy(gamma->tmpMean, x);
  int K = gamma->beta->size;
  double ans = 0.0;
  for(int k=0; k < K; k++) {
    double xy = exp(mcmclib_mvnorm_lpdf_compute(gamma->qdk[k], y));
    ans +=  xy * gsl_vector_get(gamma->weights, k);
  }
  return log(ans);
}

mcmclib_at7_suff* mcmclib_at7_suff_alloc(mcmclib_at7_gamma* g, int t0) {
  return (mcmclib_at7_suff*) mcmclib_mixem_online_alloc(g->mu, g->Sigma, g->beta, 0.5, t0);
}

void mcmclib_at7_suff_free(void* p) {
  mcmclib_mixem_online_free((mcmclib_mixem_online*) p);
}

mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p logdistr, void* logdistr_data,
			       gsl_vector* x, int t0,
			       const gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat){
  mcmclib_at7_gamma* gamma = mcmclib_at7_gamma_alloc(beta_hat,
						     mu_hat,
						     Sigma_hat);
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, at7_q_sample, at7_q_d,
				       gamma, mcmclib_at7_gamma_free);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  mcmclib_at7_suff* suff = mcmclib_at7_suff_alloc(gamma, t0);
  mcmclib_amh* ans = mcmclib_amh_alloc(mh, suff, mcmclib_at7_update,
				       mcmclib_at7_suff_free);
  return ans;
}

void mcmclib_at7_set_sf_all(mcmclib_amh* p, double sf) {
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
	  gsl_matrix_get(AT7_GAMMA(p)->Sigma[1], 0, 0))			\

void mcmclib_at7_update(mcmclib_amh* p) {
  mcmclib_mixem_online* em = (mcmclib_mixem_online*) AT7_SUFF(p);
  mcmclib_mixem_online_update(em, p->mh->x);
  if((p->n) > em->n0) {
    at7_gamma_update_Sigma(AT7_GAMMA(p));
  }
}
