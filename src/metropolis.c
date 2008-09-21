#include "metropolis.h"

int mcmclib_metropolis_symmetric_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x, distrfun_p logdistr, void* data) {
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = logdistr(data, old);
	if(!isfinite(loglik_old))
		return 1;

	loglik_new = logdistr(data, x);
	if(isfinite(loglik_new) && (loglik_new >= loglik_old))
		return 1;

	lik_ratio = exp(loglik_new - loglik_old);
	if(isfinite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio))
		return 1;

	gsl_vector_memcpy(x, old);
	return 0;
}
