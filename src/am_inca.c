#include "am_inca.h"
#include "gauss_am.h"

mcmclib_am_inca* mcmclib_am_inca_alloc(gsl_rng* r,
				       distrfun_p logdistr, void* data,
				       gsl_vector** x, int M, int t0,
				       const gsl_matrix* sigma_prop) {
  mcmclib_am_inca* a = (mcmclib_am_inca*) malloc(sizeof(mcmclib_am_inca));
  a->inca = mcmclib_inca_alloc((mcmclib_amh*) mcmclib_gauss_am_alloc(r, logdistr, data,
								     gsl_vector_alloc(1),
								     sigma_prop, t0 * M),
			       x, M);
  return a;
}

void mcmclib_am_inca_free(mcmclib_am_inca* p) {
  gsl_vector_free(p->inca->amh->mh->x);
  mcmclib_gauss_am_free((mcmclib_gauss_am*) p->inca->amh);
  free(p);
}

int mcmclib_am_inca_update(mcmclib_am_inca* p) {
  return mcmclib_inca_update(p->inca);
}
