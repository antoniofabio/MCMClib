#include "gauss_am.h"
#include "vector_stats.h"
#include "metropolis.h"

mcmclib_gauss_am* mcmclib_gauss_am_alloc(const gsl_matrix* sigma_zero, int t0) {
	mcmclib_gauss_am* ans = (mcmclib_gauss_am*) malloc(sizeof(mcmclib_gauss_am));
	int d = sigma_zero->size1;
	ans->sigma_zero = gsl_matrix_alloc(d, d);
	ans->t0 = t0;
	ans->mean = gsl_vector_alloc(d);
	ans->cov = gsl_matrix_alloc(d, d);
	ans->t = 0;
	ans->old = gsl_vector_alloc(d);
	gsl_matrix_memcpy(ans->sigma_zero, sigma_zero);
	ans->sf = (2.38 * 2.38) / d;
	ans->sigma_prop = gsl_matrix_alloc(d, d);

	gsl_vector_set_zero(ans->mean);
	gsl_matrix_set_zero(ans->cov);
	return ans;
}

void mcmclib_gauss_am_free(mcmclib_gauss_am* p) {
	gsl_matrix_free(p->sigma_prop);
	gsl_matrix_free(p->sigma_zero);
	gsl_vector_free(p->mean);
	gsl_matrix_free(p->cov);
	gsl_vector_free(p->old);
	free(p);
}

int mcmclib_gauss_am_update(mcmclib_gauss_am* e, const gsl_rng* r,
	distrfun_p logdistr, gsl_vector* x, void* data) {

	int d = x->size;
	gsl_vector* old = e->old;
	gsl_matrix* sigma_zero = e->sigma_zero;
	gsl_matrix* cov = e->cov;
	gsl_matrix* sigma_prop =  e->sigma_prop;
	gsl_vector* mean = e->mean;
	int t0 = e->t0;
	int *t = &(e->t);

	gsl_vector_memcpy(old, x);
	/*sample a new value for 'x', with mean 'old' and assigned covariance matrix*/
	mcmclib_mvnorm(r, ((*t)+1) < t0 ? sigma_zero : sigma_prop, x);
	gsl_vector_add(x, old);

	mcmclib_metropolis_symmetric_step(r, old, x, logdistr, data);

	/*adapt extra parameters*/
	mcmclib_covariance_update(cov, mean, t, x);
	gsl_matrix_memcpy(sigma_prop, cov);
	gsl_matrix_scale(sigma_prop, e->sf);

	return 0;
}
