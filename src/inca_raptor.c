#include "inca_raptor.h"
#include "raptor.h"

mcmclib_inca* mcmclib_inca_raptor_alloc(gsl_rng* r,
					distrfun_p logdistr, void* logdistr_data,
					gsl_vector** x, int M,
					int t0, gsl_matrix* Sigma_zero,
					gsl_vector* beta_hat,
					gsl_vector** mu_hat,
					gsl_matrix** Sigma_hat){
  gsl_vector* x_wrk = gsl_vector_alloc(x[0]->size);
  mcmclib_amh* amh = mcmclib_raptor_alloc(r, logdistr, logdistr_data,
					  x_wrk, t0, Sigma_zero,
					  beta_hat, mu_hat, Sigma_hat);
  return mcmclib_inca_alloc(amh, x, M);
}

void mcmclib_inca_raptor_free(mcmclib_inca* p) {
  gsl_vector_free(p->amh->mh->x);
  mcmclib_raptor_free(p->amh);
  mcmclib_inca_free(p);
}
