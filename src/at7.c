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

#define AT7_GAMMA(p) ((mcmclib_at7_gamma*)((mcmclib_amh*) p)->mh->q->gamma)
#define AT7_SUFF(p) ((mcmclib_at7_suff*)((mcmclib_amh*) p)->suff)

mcmclib_at7_gamma* mcmclib_at7_gamma_alloc(const gsl_vector* beta,
					   gsl_vector** mu,
					   gsl_matrix** Sigma) {
  mcmclib_at7_gamma* ans = (mcmclib_at7_gamma*) malloc(sizeof(mcmclib_at7_gamma));
  gsl_vector_memcpy(ans->beta, beta);
  for(int k=0; k < beta->size; k++) {
    gsl_vector_memcpy(ans->mu[k], mu[k]);
    gsl_matrix_memcpy(ans->Sigma[k], Sigma[k]);
    ans->pik[k] = mcmclib_mvnorm_lpdf_alloc(ans->mu[k], ans->Sigma[k]->data);
  }
  gsl_vector_memcpy(ans->beta, beta);
  ans->pi = mcmclib_mixnorm_lpdf_alloc(ans->beta, ans->pik);
  return ans;
}

void mcmclib_at7_gamma_free(mcmclib_at7_gamma* p) {
  mcmclib_mixnorm_lpdf_free(p->pi);
  gsl_vector_free(p->beta);
  for(int k=0; k < p->beta->size; k++) {
    mcmclib_mvnorm_lpdf_free(p->pik[k]);
    gsl_vector_free(p->mu[k]);
    gsl_matrix_free(p->Sigma[k]);
  }
  free(p);
}

static void at7_components_weights(mcmclib_at7_gamma* g, gsl_vector* x) {
  mcmclib_mvnorm_lpdf** pik = g->pik;
  double sum=0.0;
  for(int k=0; k < x->size; k++) {
    double nuk = exp(mcmclib_mvnorm_lpdf_compute(pik[k], x));
    gsl_vector_set(x, k, nuk);
    sum += nuk;
  }
  gsl_vector_scale(x, 1.0 / sum);
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
  at7_components_weights(gamma, gamma->weights);
  int k = sample(q->r, gamma->weights);
  gsl_vector* oldx = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(oldx, x);
  mcmclib_mvnorm(q->r, gamma->Sigma[k], x);
  gsl_vector_add(x, oldx);
  gsl_vector_free(oldx);
}

static double at7_q_d(void* in_gamma, gsl_vector* x, gsl_vector* y) {
  mcmclib_at7_gamma* gamma = (mcmclib_at7_gamma*) in_gamma;
  for(int k=0; k < gamma->beta->size; k++)
    gsl_vector_memcpy(gamma->mu[k], x);
  return mcmclib_mixnorm_lpdf_compute(gamma->pi, y);
}

static mcmclib_mh_q* at7_q_alloc(gsl_rng* r, distrfun_p logdistr, void* logdistr_data,
				 gsl_matrix* Sigma0, mcmclib_at7_gamma* gamma) {
  return mcmclib_mh_q_alloc(r, at7_q_sample, gamma,
			    at7_q_d, gamma,
			    gamma);
}

mcmclib_at7_suff* mcmclib_at7_suff_alloc(mcmclib_at7_gamma* g,
					 const gsl_vector* beta_hat,
					 gsl_vector** mu_hat,
					 gsl_matrix** Sigma_hat,
					 int t0) {
  mcmclib_at7_suff* a = (mcmclib_at7_suff*) malloc(sizeof(mcmclib_at7_suff));
  int K = g->beta->size;
  int dim = g->mu[0]->size;
  a->beta_hat = gsl_vector_alloc(K);
  gsl_vector_memcpy(a->beta_hat, beta_hat);
  a->mu_hat = malloc(K * sizeof(gsl_vector*));
  a->Sigma_hat = malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    a->mu_hat[k] = gsl_vector_alloc(dim);
    gsl_vector_memcpy(a->mu_hat[k], mu_hat[k]);
    a->Sigma_hat[k] = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(a->Sigma_hat[k], Sigma_hat[k]);
  }
  a->em = mcmclib_mixem_online_alloc(a->mu_hat, a->Sigma_hat, a->beta_hat, 0.5, t0);
  return a;
}

void mcmclib_at7_suff_free(mcmclib_at7_suff* p) {
  mcmclib_mixem_online_free(p->em);
  for(int k=0; k < p->beta_hat->size; k++) {
    gsl_vector_free(p->mu_hat[k]);
    gsl_matrix_free(p->Sigma_hat[k]);
  }
  free(p->Sigma_hat);
  free(p->mu_hat);
  gsl_vector_free(p->beta_hat);
  free(p);
}

mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p logdistr, void* logdistr_data,
			       gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
			       const gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat){
  mcmclib_at7_gamma* gamma = mcmclib_at7_gamma_alloc(beta_hat,
						     mu_hat,
						     Sigma_hat);
  mcmclib_mh_q* q = at7_q_alloc(r, logdistr, logdistr_data,
				Sigma_zero, gamma);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  mcmclib_at7_suff* suff = mcmclib_at7_suff_alloc(gamma,
						  beta_hat,
						  mu_hat,
						  Sigma_hat,
						  t0);
  mcmclib_amh* ans = mcmclib_amh_alloc(mh, suff, mcmclib_at7_update);

  mcmclib_at7_set_sf_all(ans, 2.38*2.38 / (double) x->size);
  return ans;
}

void mcmclib_at7_free(mcmclib_amh* p) {
  mcmclib_at7_suff_free(AT7_SUFF(p));
  mcmclib_at7_gamma_free(AT7_GAMMA(p));
  //mcmclib_at7_q_free((mcmclib_at7_q*) p->mh->q); FIXME
  mcmclib_mh_free(p->mh);
  mcmclib_amh_free(p);
}

void mcmclib_at7_set_sf_all(mcmclib_amh* p, double sf) {
  gsl_vector_set_all(AT7_SUFF(p)->scaling_factors, sf);
}

void mcmclib_at7_set_sf(mcmclib_amh* p, const gsl_vector* sf) {
  gsl_vector_memcpy(AT7_SUFF(p)->scaling_factors, sf);
}

void mcmclib_at7_update(void* in_p) {
  mcmclib_amh* p = (mcmclib_amh*) in_p;
  mcmclib_at7_suff* s = AT7_SUFF(p);
  mcmclib_mixem_online* em = s->em;
  mcmclib_mixem_online_update(em, p->mh->x);
  if((p->n) <= em->n0)
    return;
  //TODO
}
