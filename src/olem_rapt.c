#include <gsl/gsl_math.h>
#include "olem_rapt.h"
#include "region_mixnorm.h"
#include "mvnorm.h"

#define MOR mcmclib_olemrapt

MOR* mcmclib_olemrapt_alloc(gsl_rng* r,
			    distrfun_p logdistr, void* logdistr_data,
			    gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
			    gsl_vector* beta_hat,
			    gsl_vector** mu_hat,
			    gsl_matrix** Sigma_hat){
  int dim = x->size;
  MOR* a = (MOR*) malloc(sizeof(MOR));
  int K = beta_hat->size;
  a->beta_hat = beta_hat;
  a->mu_hat = mu_hat;
  a->Sigma_hat = Sigma_hat;

  a->pik_hat = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++)
    a->pik_hat[k] = mcmclib_mvnorm_lpdf_alloc(mu_hat[k], Sigma_hat[k]->data);
  a->pi_hat = mcmclib_mixnorm_lpdf_alloc(a->beta_hat, a->pik_hat);

  gsl_matrix** tmp = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    tmp[k] = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(tmp[k], Sigma_hat[k]);
    gsl_matrix_scale(tmp[k], 2.38 * 2.38 / ((double) dim));
  }
  a->rapt = mcmclib_rapt_alloc(r, logdistr, logdistr_data,
			       x, t0, Sigma_zero,
			       K, tmp,
			       mcmclib_region_mixnorm_compute, a->pi_hat);
  for(int k=0; k<K; k++)
    gsl_matrix_free(tmp[k]);
  free(tmp);

  a->em = mcmclib_mixem_online_alloc(mu_hat, Sigma_hat, beta_hat, 0.5, t0);

  return a;
}

void mcmclib_olemrapt_free(MOR* p) {
  mcmclib_mixem_online_free(p->em);
  mcmclib_rapt_free(p->rapt);
  mcmclib_mixnorm_lpdf_free(p->pi_hat);
  for(int k=0; k< p->beta_hat->size; k++)
    mcmclib_mvnorm_lpdf_free(p->pik_hat[k]);
  free(p->pik_hat);
  free(p);
}

int mcmclib_olemrapt_update(MOR* p) {
  int ans = mcmclib_rapt_update(p->rapt);
  mcmclib_mixem_online_update(p->em, p->rapt->current_x);
  return ans;
}

void mcmclib_olemrapt_update_proposals(MOR* p) {
  mcmclib_rapt* r = p->rapt;
  if((r->t) <= r->t0)
    return;
  mcmclib_rapt_update_proposals_custom(r, p->Sigma_hat, r->global_variance);
}
