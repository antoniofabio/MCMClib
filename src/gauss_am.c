#include "gauss_am.h"
#include "vector_stats.h"

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

	gsl_vector_set_zero(ans->mean);
	gsl_matrix_set_zero(ans->cov);
	return ans;
}

void mcmclib_gauss_am_free(mcmclib_gauss_am* p) {
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
	gsl_vector* mean = e->mean;
	int t0 = e->t0;
	int *t = &(e->t);

	(*t)++;
	gsl_vector_memcpy(old, x);
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = logdistr(data, x);

	/*sample a new value for 'x', with mean 'old' and assigned covariance matrix*/
	mcmclib_mvnorm(r, (*t) < t0 ? sigma_zero : cov, x);
	gsl_vector_add(x, old);

	if(!isfinite(loglik_old))
		return 0;

	loglik_new = logdistr(data, x);
	if(loglik_new >= loglik_old)
		return 0;

	lik_ratio = exp(loglik_new - loglik_old);
	if(isfinite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio))
		return 0;

	gsl_vector_memcpy(x, old);

	/*adapt extra parameters*/
	(*t)--;
	mcmclib_covariance_update(cov, mean, t, x);

	return 0;
}
