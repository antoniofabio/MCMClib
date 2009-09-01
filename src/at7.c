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

mcmclib_at7_gamma* mcmclib_at7_gamma_alloc(gsl_vector* beta_hat,
					   gsl_vector** mu_hat,
					   gsl_matrix** Sigma_hat) {
  return (mcmclib_at7_gamma*) mcmclib_raptor_gamma_alloc(beta_hat, mu_hat, Sigma_hat);
}

void mcmclib_at7_gamma_free(mcmclib_at7_gamma* p) {
  mcmclib_raptor_gamma_free((mcmclib_raptor_gamma*) p);
}

mcmclib_at7_suff* mcmclib_at7_suff_alloc(mcmclib_at7_gamma* g, int t0) {
  mcmclib_at7_suff* a = (mcmclib_at7_suff*) malloc(sizeof(mcmclib_at7_suff));
  return a;
}

void mcmclib_at7_suff_free(mcmclib_at7_suff* p) {
  free(p);
}

static void at7_regions_weights(mcmclib_at7_gamma* g, gsl_vector* x) {
  mcmclib_mvnorm_lpdf** pik_hat = g->pik_hat;
  for(int k=0; k < x->size; k++)
    gsl_vector_set(x, k, exp(mcmclib_mvnorm_lpdf_compute(pik_hat[k], x)));
}

mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p logdistr, void* logdistr_data,
			       gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
			       gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat){
  //int K = beta_hat->size;

  mcmclib_at7_gamma* gamma = mcmclib_at7_gamma_alloc(beta_hat,
						     mu_hat,
						     Sigma_hat);
  /*  mcmclib_mh_q* q = mcmclib_at7_q_alloc(r, logdistr, logdistr_data,
					Sigma_zero, K, Sigma_hat,
					gamma);*/ //FIXME
  mcmclib_mh_q* q = NULL;
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  mcmclib_at7_suff* suff = mcmclib_at7_suff_alloc(gamma, t0);
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
  gsl_vector_set_all(AT7_GAMMA(p)->scaling_factors, sf);
}

void mcmclib_at7_set_sf(mcmclib_amh* p, const gsl_vector* sf) {
  gsl_vector_memcpy(AT7_GAMMA(p)->scaling_factors, sf);
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
