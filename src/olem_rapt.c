#include <gsl/gsl_math.h>
#include "olem_rapt.h"
#include "region_mixnorm.h"
#include "mvnorm.h"

#define MOR mcmclib_olemrapt

MOR* mcmclib_olemrapt_alloc(gsl_rng* r,
	       distrfun_p logdistr, void* logdistr_data,
	       gsl_vector* x, int t0,
	       gsl_vector* beta_hat,
	       gsl_vector** mu_hat,
	       gsl_matrix** Sigma_hat){
  int dim = x->size;
  MOR* a = (MOR*) malloc(sizeof(MOR));
  int K = beta_hat->size;

  gsl_matrix* Sigma_whole = gsl_matrix_alloc(dim, dim);
  /**FIXME: compute global variance from group means and variances*/
  gsl_matrix_set_identity(Sigma_whole);

  a->pik_hat = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++)
    a->pik_hat[k] = mcmclib_mvnorm_lpdf_alloc(mu_hat[k], Sigma_hat[k]->data);
  a->pi_hat = mcmclib_mixnorm_lpdf_alloc(a->beta_hat, a->pik_hat);

  a->rapt = mcmclib_rapt_alloc(r, logdistr, logdistr_data,
			       x, t0,
			       Sigma_whole,
			       K,
			       Sigma_hat,
			       mcmclib_region_mixnorm_compute, a->pi_hat);
  gsl_matrix_free(Sigma_whole);

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
  return mcmclib_rapt_update(p->rapt);
}

void mcmclib_olemrapt_update_proposals(MOR* p) {
  mcmclib_rapt* r = p->rapt;
  if((r->t) <= r->t0)
    return;
  double sf = 2.38 * 2.38 / ((double) p->mu_hat[0]->size);
  for(int k=0; k < r->K; k++) {
    gsl_matrix_memcpy(r->sigma_local[k], p->Sigma_hat[k]);
    gsl_matrix_add(r->sigma_local[k], r->Sigma_eps);
    gsl_matrix_scale(r->sigma_local[k], sf);
  }
  gsl_matrix_memcpy(r->sigma_whole, r->global_variance);
  gsl_matrix_add(r->sigma_whole, r->Sigma_eps);
  gsl_matrix_scale(r->sigma_whole, sf);
}
