#include "gauss_rw.h"

mcmclib_gauss_rw_data* mcmclib_gauss_rw_alloc(double step_size, int dim) {
	mcmclib_gauss_rw_data* ans = (mcmclib_gauss_rw_data*) malloc(sizeof(mcmclib_gauss_rw_data));
	ans->old = gsl_vector_alloc(dim);
	ans->step_size = step_size;
	return ans;
}

void mcmclib_gauss_rw_free(mcmclib_gauss_rw_data* p) {
	gsl_vector_free(p->old);
	free(p);
}

int mcmclib_gauss_rw(mcmclib_gauss_rw_data* e, const gsl_rng* r,
	distrfun_p logdistr, gsl_vector* x, void* data) {
	int n = e->old->size;
	double step_size = e->step_size;
	gsl_vector* old = e->old;
	gsl_vector_memcpy(old, x);
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = logdistr(data, x);

	while(n--) {
		gsl_vector_set(x, n,
			gsl_vector_get(x, n) + gsl_ran_gaussian(r, step_size));
	}

	if(!isfinite(loglik_old))
		return 0;

	loglik_new = logdistr(data, x);
	if(loglik_new >= loglik_old)
		return 0;

	lik_ratio = exp(loglik_new - loglik_old);
	if(isfinite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio))
		return 0;

	gsl_vector_memcpy(x, old);
	return 0;
}
