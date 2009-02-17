#include "inca.h"

mcmclib_inca* mcmclib_inca_alloc(gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 proposal_distr_t qd_fun, void* qd_data,
				 samplerfun_p q_sampler, void* q_sampler_data,
				 mcmclib_amh_update_gamma_p update_gamma,
				 void* update_gamma_data,
				 gsl_vector** x, int M) {
  mcmclib_inca* p = (mcmclib_inca*) malloc(sizeof(mcmclib_inca));
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data,
				    gsl_vector_alloc(x[0]->size),
				    qd_fun, qd_data,
				    q_sampler, q_sampler_data);
  p->amh = mcmclib_amh_alloc(mh, update_gamma, update_gamma_data);
  p->x = x;
  p->M = M;
  return p;
}

void mcmclib_inca_free(mcmclib_inca* p) {
  gsl_vector_free(p->amh->mh->x);
  mcmclib_mh_free(p->amh->mh);
  mcmclib_amh_free(p->amh);
  free(p);
}

int mcmclib_inca_update(mcmclib_inca* p) {
  mcmclib_amh* amh = p->amh;
  mcmclib_mh* mh = amh->mh;
  gsl_vector** x = p->x;
  for(int m=0; m < p->M; m++) {
    gsl_vector_memcpy(mh->x, x[m]);
    mcmclib_amh_update(amh);
    gsl_vector_memcpy(x[m], mh->x);
  }
  return(mh->last_accepted);
}
