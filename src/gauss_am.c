#include "gauss_am.h"
#include "vector_stats.h"

mcmclib_gauss_am_data* mcmclib_gauss_am_alloc(const gsl_matrix* sigma_zero, int t0) {
	mcmclib_gauss_am_data* ans = (mcmclib_gauss_am_data*) malloc(sizeof(mcmclib_gauss_am_data));
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

void mcmclib_gauss_am_free(mcmclib_gauss_am_data* p) {
	gsl_matrix_free(p->sigma_zero);
	gsl_vector_free(p->mean);
	gsl_matrix_free(p->cov);
	gsl_vector_free(p->old);
	free(p);
}

int mcmclib_gauss_am(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	mcmclib_gauss_am_data* e) {

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

	loglik_old = loglik(x, data);

	/*sample a new value for 'x', with mean 'old' and assigned covariance matrix*/
	mcmclib_mvnorm(r, (*t) < t0 ? sigma_zero : cov, x);
	gsl_vector_add(x, old);

	if(!isfinite(loglik_old)) {
		gsl_vector_free(old);
		return 0;
	}

	loglik_new = loglik(x, data);
	if(loglik_new >= loglik_old) {
		gsl_vector_free(old);
		return 0;
	}

	lik_ratio = exp(loglik_new - loglik_old);
	if(isfinite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio)) {
		gsl_vector_free(old);
		return 0;
	}

	gsl_vector_memcpy(x, old);
	gsl_vector_free(old);

	/*adapt extra parameters*/
	(*t)--;
	mcmclib_covariance_update(cov, mean, t, x);

	return 0;
}
