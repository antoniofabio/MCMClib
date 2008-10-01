#include "gauss_am.h"
#include "vector_stats.h"
#include "metropolis.h"

mcmclib_gauss_am* mcmclib_gauss_am_alloc(gsl_rng* r,
	distrfun_p logdistr, void* logdistr_data, gsl_vector* start_x,
	const gsl_matrix* sigma_zero, int t0) {
	mcmclib_gauss_am* ans = (mcmclib_gauss_am*) malloc(sizeof(mcmclib_gauss_am));
	int d = sigma_zero->size1;
	ans->r = r;
	ans->logdistr = logdistr;
	ans->logdistr_data = logdistr_data;
	ans->current_x = gsl_vector_alloc(d);
	gsl_vector_memcpy(ans->current_x, start_x);
	ans->old = gsl_vector_alloc(d);

	ans->sigma_zero = gsl_matrix_alloc(d, d);
	ans->t0 = t0;
	ans->mean = gsl_vector_alloc(d);
	ans->cov = gsl_matrix_alloc(d, d);
	ans->t = 0;
	gsl_matrix_memcpy(ans->sigma_zero, sigma_zero);
	ans->sf = (2.38 * 2.38) / d;
	ans->sigma_prop = gsl_matrix_alloc(d, d);

	gsl_vector_set_zero(ans->mean);
	gsl_matrix_set_zero(ans->cov);
	return ans;
}

void mcmclib_gauss_am_free(mcmclib_gauss_am* p) {
	gsl_vector_free(p->current_x);
	gsl_vector_free(p->old);
	gsl_matrix_free(p->sigma_prop);
	gsl_matrix_free(p->sigma_zero);
	gsl_vector_free(p->mean);
	gsl_matrix_free(p->cov);
	free(p);
}

int mcmclib_gauss_am_update(mcmclib_gauss_am* p) {
	gsl_rng* r = p->r;
	distrfun_p logdistr = p->logdistr;
	void* data = p->logdistr_data;
	gsl_vector* x = p->current_x;
	int d = x->size;
	gsl_vector* old = p->old;
	gsl_matrix* sigma_zero = p->sigma_zero;
	gsl_matrix* cov = p->cov;
	gsl_matrix* sigma_prop =  p->sigma_prop;
	gsl_vector* mean = p->mean;
	int t0 = p->t0;
	int *t = &(p->t);

	gsl_vector_memcpy(old, x);
	/*sample a new value for 'x', with mean 'old' and assigned covariance matrix*/
	mcmclib_mvnorm(r, ((*t)+1) < t0 ? sigma_zero : sigma_prop, x);
	gsl_vector_add(x, old);

	mcmclib_metropolis_symmetric_step(r, old, x, logdistr, data);

	/*adapt extra parameters*/
	mcmclib_covariance_update(cov, mean, t, x);
	gsl_matrix_memcpy(sigma_prop, cov);
	gsl_matrix_scale(sigma_prop, p->sf);

	return 0;
}
