#include "gauss_inca.h"
#include "vector_stats.h"

/*FIXME*/
mcmclib_gauss_inca* mcmclib_gauss_inca_alloc(mcmclib_gauss_inca_pool* pool) {
	return NULL;
}

void mcmclib_gauss_inca_free(mcmclib_gauss_inca* p) {
	if(p)
		free(p);
}

/*FIXME*/
int mcmclib_gauss_inca_update(mcmclib_gauss_inca* e, const gsl_rng* r,
	distrfun_p logdistr, gsl_vector* x, void* data) {
	return 0;
}
