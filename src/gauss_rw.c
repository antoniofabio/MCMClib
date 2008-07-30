#include "gauss_rw.h"

int mcmclib_gauss_rw(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	const double step_size) {
	int n = x->size;
	gsl_vector* old = gsl_vector_alloc(n);
	gsl_vector_memcpy(old, x);
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = loglik(x, data);

	while(n--) {
		gsl_vector_set(x, n,
			gsl_vector_get(x, n) + gsl_ran_gaussian(r, step_size));
	}

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
	return 0;
}
