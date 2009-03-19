#include "inca_rapt.h"
#include "rapt.h"

mcmclib_inca* mcmclib_inca_rapt_alloc(gsl_rng* r,
				      distrfun_p logdistr, void* logdistr_data,
				      gsl_vector** x, int M, int t0,
				      const gsl_matrix* sigma_whole, int K,
				      gsl_matrix** sigma_local,
				      region_fun_t which_region,
				      void* which_region_data){
  gsl_vector* x_wrk = gsl_vector_alloc(x[0]->size);
  mcmclib_amh* amh = mcmclib_rapt_alloc(r, logdistr, logdistr_data, x_wrk,
					t0, sigma_whole, K, sigma_local,
					which_region, which_region_data);
  return mcmclib_inca_alloc(amh, x, M);
}

void mcmclib_inca_rapt_free(mcmclib_inca* p) {
  gsl_vector_free(p->amh->mh->x);
  mcmclib_rapt_free(p->amh);
  mcmclib_inca_free(p);
}
