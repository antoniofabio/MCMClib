#include <gsl/gsl_math.h>
#include "raptor.h"
#include "region_mixnorm.h"
#include "mvnorm.h"

mcmclib_raptor_gamma* mcmclib_raptor_gamma_alloc(gsl_vector* beta_hat,
						 gsl_vector** mu_hat,
						 gsl_matrix** Sigma_hat) {
  mcmclib_raptor_gamma* a = (mcmclib_raptor_gamma*) malloc(sizeof(mcmclib_raptor_gamma));
  a->beta_hat = gsl_vector_alloc(beta_hat->size);
  gsl_vector_memcpy(a->beta_hat, beta_hat);
  int K = beta_hat->size;
  int d = mu_hat[0]->size;
  a->mu_hat = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  a->Sigma_hat = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k < K; k++) {
    a->mu_hat[k] = gsl_vector_alloc(d);
    gsl_vector_memcpy(a->mu_hat[k], mu_hat[k]);
    a->Sigma_hat[k] = gsl_matrix_alloc(d, d);
    gsl_matrix_memcpy(a->Sigma_hat[k], Sigma_hat[k]);
  }

  a->pik_hat = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++)
    a->pik_hat[k] = mcmclib_mvnorm_lpdf_alloc(a->mu_hat[k], a->Sigma_hat[k]->data);

  a->pi_hat = mcmclib_mixnorm_lpdf_alloc(a->beta_hat, a->pik_hat);
  return a;
}

void mcmclib_raptor_gamma_free(mcmclib_raptor_gamma* p) {
  mcmclib_mixnorm_lpdf_free(p->pi_hat);
  for(int k=0; k< p->beta_hat->size; k++) {
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

mcmclib_raptor_suff* mcmclib_raptor_suff_alloc(mcmclib_raptor_gamma* g, int t0) {  
  return (mcmclib_raptor_suff*) mcmclib_mixem_online_alloc(g->mu_hat, g->Sigma_hat, g->beta_hat, 0.5, t0);
}

void mcmclib_raptor_suff_free(mcmclib_raptor_suff* p) {
  mcmclib_mixem_online_free((mcmclib_mixem_online*) p);
}

int raptor_which_region_fun(gsl_vector* x, void* in_g) {
  mcmclib_raptor_gamma* g = (mcmclib_raptor_gamma*) in_g;
  return mcmclib_region_mixnorm_compute(x, g->pi_hat);
}

void mcmclib_raptor_update(void* in_p);

mcmclib_amh* mcmclib_raptor_alloc(gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
				 gsl_vector* beta_hat,
				 gsl_vector** mu_hat,
				 gsl_matrix** Sigma_hat){
  int K = beta_hat->size;

  mcmclib_raptor_gamma* gamma = mcmclib_raptor_gamma_alloc(beta_hat,
							   mu_hat,
							   Sigma_hat);
  mcmclib_mh_q* q = mcmclib_rapt_q_alloc(r, logdistr, logdistr_data,
					 Sigma_zero, K, Sigma_hat,
					 raptor_which_region_fun, gamma);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  mcmclib_raptor_suff* suff = mcmclib_raptor_suff_alloc(gamma, t0);

  return mcmclib_amh_alloc(mh, suff, mcmclib_raptor_update);
}

void mcmclib_raptor_free(mcmclib_amh* p) {
  mcmclib_raptor_suff_free((mcmclib_raptor_suff*) p->suff);
  mcmclib_raptor_gamma_free((mcmclib_raptor_gamma*) p->mh->q->gamma);
  mcmclib_rapt_q_free(p->mh->q);
  mcmclib_mh_free(p->mh);
  mcmclib_amh_free(p);
}

void mcmclib_raptor_update(void* in_p) {
  mcmclib_amh* p = (mcmclib_amh*) in_p;
  mcmclib_mixem_online* s = (mcmclib_mixem_online*) p->suff;
  mcmclib_mixem_online_update(s, p->mh->x);
  if((p->n) <= s->n0)
    return;
  mcmclib_rapt_update_proposals_custom(p, s->Sigma, s->Sigma_global);
}
