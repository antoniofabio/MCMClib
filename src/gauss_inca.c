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
	ans->id = 0;
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

mcmclib_gauss_inca* mcmclib_gauss_inca_alloc(mcmclib_gauss_inca_pool* p) {
	if(p->id == p->K)
		return NULL;
	mcmclib_gauss_inca* ans = (mcmclib_gauss_inca*) malloc(sizeof(mcmclib_gauss_inca));
	int d = (p->Sigma_zero)->size1;
	ans->old = gsl_vector_alloc(d);
	ans->id = p->id;
	(p->id)++;
	return ans;
}

void mcmclib_gauss_inca_free(mcmclib_gauss_inca* p) {
	if(p) {
		gsl_vector_free(p->old);
		free(p);
	}
}

/**returns 1 for accept, 0 for reject
@param r GSL RNG
@param old old value
@param x vector holding current value (will be eventually updated!)
@param logdistr ptr to log-distribution function
@param data extra data for 'logdistr'
*/
int metropolis_symmetric_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x, distrfun_p logdistr, void* data) {
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = logdistr(data, old);
	if(!isfinite(loglik_old))
		return 1;

	loglik_new = logdistr(data, x);
	if(loglik_new >= loglik_old)
		return 1;

	lik_ratio = exp(loglik_new - loglik_old);
	if(isfinite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio))
		return 1;

	gsl_vector_memcpy(x, old);
	return 0;
}

/*FIXME: currently just ignoring other chains*/
int mcmclib_gauss_inca_update(mcmclib_gauss_inca* e, const gsl_rng* r,
	distrfun_p logdistr, gsl_vector* x, void* data) {

	int d = x->size;
	gsl_vector* old = e->old;
	gsl_matrix* sigma_zero = e->p->Sigma_zero;
	int id = e->id;
	gsl_matrix* cov = (e->p->variance)[id];
	gsl_vector* mean = (e->p->mean)[id];
	int t0 = e->p->t0;
	int *t = (e->p->t) + id;

	gsl_vector_memcpy(old, x);
	mcmclib_mvnorm(r, ((*t)+1) < t0 ? sigma_zero : cov, x);
	gsl_vector_add(x, old);

	metropolis_symmetric_step(r, old, x, logdistr, data);

	/*adapt extra parameters*/
	mcmclib_covariance_update(cov, mean, t, x);

	return 0;
}
