#include "olemrapt_inca.h"

#define MOI mcmclib_olemrapt_inca

MOI* mcmclib_olemrapt_inca_alloc(gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 gsl_vector** x, int t0, gsl_matrix* Sigma_zero,
				 gsl_vector* beta_hat,
				 gsl_vector** mu_hat,
				 gsl_matrix** Sigma_hat,
				 int M) {
  MOI* a = (MOI*) malloc(sizeof(MOI));
  a->M = M;
  a->ss = (mcmclib_olemrapt**) malloc(M * sizeof(mcmclib_olemrapt*));
  for(int m=0; m<M; m++)
    a->ss[m] = mcmclib_olemrapt_alloc(r, logdistr, logdistr_data,
				      x[m], t0, Sigma_zero,
				      beta_hat, mu_hat, Sigma_hat);
  a->s = mcmclib_mixolem_suff_alloc(beta_hat->size, mu_hat[0]->size);
  a->workspace = mcmclib_mixolem_suff_alloc(beta_hat->size, mu_hat[0]->size);
  a->gamma_hat = mcmclib_mixolem_suff_alloc(beta_hat->size, mu_hat[0]->size);
  return a;
}

void mcmclib_olemrapt_inca_free(MOI* p){
  mcmclib_mixolem_suff_free(p->gamma_hat);
  mcmclib_mixolem_suff_free(p->workspace);
  mcmclib_mixolem_suff_free(p->s);
  for(int m=0; m < p->M; m++)
    mcmclib_olemrapt_free(p->ss[m]);
  free(p->ss);
  free(p);
}

int mcmclib_olemrapt_inca_update(mcmclib_olemrapt_inca* p) {
  int M = p->M;
  for(int m=0; m < M; m++) {
    mcmclib_olemrapt_update(p->ss[m]);
    mcmclib_mixolem_suff_memcpy(p->ss[(m+1) % M]->em->s, p->ss[m]->em->s);
    mcmclib_mixolem_suff_memcpy(p->ss[(m+1) % M]->em->gamma, p->ss[m]->em->gamma);
  }
  return 0;
}

void mcmclib_olemrapt_inca_update_proposals(mcmclib_olemrapt_inca* p){
  int M = p->M;

  /*update chains proposals*/
  for(int m=0; m < M; m++) {
    mcmclib_olemrapt_update_proposals(p->ss[m]);
    if(m < (M-1)) {
      mcmclib_mixolem_suff_memcpy(p->ss[m+1]->em->s, p->ss[m]->em->s);
      mcmclib_mixolem_suff_memcpy(p->ss[m+1]->em->gamma, p->ss[m]->em->gamma);
    }
  }
}
