#include "olemrapt_inca.h"

/**\TODO*/
mcmclib_olemrapt_inca* mcmclib_olemrapt_inca_alloc(gsl_rng* r,
						   distrfun_p logdistr, void* logdistr_data,
						   gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
						   gsl_vector* beta_hat,
						   gsl_vector** mu_hat,
						   gsl_matrix** Sigma_hat) {
  return NULL;
}

/**\TODO*/
void mcmclib_olemrapt_inca_free(mcmclib_olemrapt_inca* p){
}

/**\TODO*/
int mcmclib_olemrapt_inca_update(mcmclib_olemrapt_inca* p) {
  return 0;
}

/**\TODO*/
void mcmclib_olemrapt_inca_update_proposals(mcmclib_olemrapt_inca* p){
}
