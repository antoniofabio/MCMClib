#include "common.h"

int mcmclib_gauss_am(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	gsl_matrix* sigma_zero, int t0,
	gsl_matrix* cov, gsl_vector* mean, int* t) {

	//FIXME

	int d = x->size;
	gsl_vector* old = gsl_vector_alloc(d);
	gsl_vector_memcpy(old, x);
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = loglik(x, data);

	//TODO: multivariate normal distribution sampler

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
	/*update mean:*/
	gsl_vector_scale(mean, (double) (*t));
	gsl_vector_add(mean, x);
	gsl_vector_scale(mean, 1.0 / (double) (*t));
	/*end update mean*/
	//TODO: recursively update sample covariance matrix 'cov'

	return 0;
}
