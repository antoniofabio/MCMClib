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
  a->em = mcmclib_mixem_online_alloc(mu_hat, Sigma_hat, beta_hat, 0.5, t0);
  a->ss = (mcmclib_olemrapt**) malloc(M * sizeof(mcmclib_olemrapt*));
  for(int m=0; m<M; m++)
    a->ss[m] = mcmclib_olemrapt_alloc(r, logdistr, logdistr_data,
				      x[m], t0, Sigma_zero,
				      beta_hat, mu_hat, Sigma_hat);
  return a;
}

void mcmclib_olemrapt_inca_free(MOI* p){
  mcmclib_mixem_online_free(p->em);
  for(int m=0; m < p->M; m++)
    mcmclib_olemrapt_free(p->ss[m]);
  free(p->ss);
  free(p);
}

int mcmclib_olemrapt_inca_update(mcmclib_olemrapt_inca* p) {
  int M = p->M;
  for(int m=0; m < M; m++) {
    mcmclib_rapt_update(p->ss[m]->rapt);
    mcmclib_mixem_online_update(p->em, p->ss[m]->rapt->current_x);
  }
  return 0;
}

/*update chains proposals*/
void mcmclib_olemrapt_inca_update_proposals(mcmclib_olemrapt_inca* p){
  for(int m=0; m < p->M; m++)
    mcmclib_olemrapt_update_proposals(p->ss[m]);
}
