#include "gauss_inca.h"
#include "vector_stats.h"

/** INCA chains shared data structure allocator
@param Sigma_zero starting covariance guess
@param t0 burn-in length
@param K number of parallel chains to account for
*/
mcmclib_gauss_inca_pool* mcmclib_gauss_inca_pool_alloc(gsl_matrix* Sigma_zero, int t0, int K) {
	mcmclib_gauss_inca_pool* ans = (mcmclib_gauss_inca_pool*) malloc(sizeof(mcmclib_gauss_inca_pool));
	int d = Sigma_zero->size1;
	ans->Sigma_zero = gsl_matrix_alloc(d, d);
	gsl_matrix_memcpy(ans->Sigma_zero, Sigma_zero);
	ans->t0 = t0;
	ans->K = K;
	ans->mean = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
	ans->variance = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
	ans->t = (int*) malloc(K * sizeof(int));
	for(int i=0; i<K; i++) {
		(ans->mean)[i] = gsl_vector_alloc(d);
		gsl_vector_set_all((ans->mean)[i], 0);
		(ans->variance)[i] = gsl_matrix_alloc(d, d);
		gsl_matrix_set_identity((ans->variance)[i]);
		(ans->t)[i] = 0;
	}
	ans->mean_global = gsl_vector_alloc(d);
	gsl_vector_set_all(ans->mean_global, 0);
	ans->variance_global = gsl_matrix_alloc(d, d);
	gsl_matrix_set_identity(ans->variance_global);
	return NULL;
}

void mcmclib_gauss_inca_pool_free(mcmclib_gauss_inca_pool* p) {
	if(p) {
		gsl_matrix_free(p->variance_global);
		gsl_vector_free(p->mean_global);
		for(int i=0; i<p->K; i++) {
			gsl_matrix_free((p->variance)[i]);
			gsl_vector_free((p->mean)[i]);
		}
		free(p->t);
		free(p->variance);
		free(p->mean);
		gsl_matrix_free(p->Sigma_zero);
		free(p);
	}
}

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
