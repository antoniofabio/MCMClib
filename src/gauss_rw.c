#include <gsl/gsl_math.h>
#include "gauss_rw.h"

mcmclib_gauss_rw* mcmclib_gauss_rw_alloc(gsl_rng* r,
	distrfun_p logdistr, void* logdistr_data, gsl_vector* start_x, double step_size) {
	int dim = start_x->size;
	mcmclib_gauss_rw* ans = (mcmclib_gauss_rw*) malloc(sizeof(mcmclib_gauss_rw));
	ans->r = r;
	ans->logdistr = logdistr;
	ans->logdistr_data = logdistr_data;
	ans->current_x = gsl_vector_alloc(dim);
	gsl_vector_memcpy(ans->current_x, start_x);
	ans->old = gsl_vector_alloc(dim);
	ans->step_size = step_size;
	return ans;
}

void mcmclib_gauss_rw_free(mcmclib_gauss_rw* p) {
	gsl_vector_free(p->current_x);
	gsl_vector_free(p->old);
	free(p);
}

int mcmclib_gauss_rw_update(mcmclib_gauss_rw* p) {
	int n = p->old->size;
	gsl_rng* r = p->r;
	double step_size = p->step_size;
	gsl_vector* old = p->old;
	gsl_vector* x = p->current_x;
	gsl_vector_memcpy(old, x);
	distrfun_p logdistr = p->logdistr;
	void* data = p->logdistr_data;
	double loglik_old, loglik_new, lik_ratio;

	loglik_old = logdistr(data, x);

	while(n--) {
		gsl_vector_set(x, n,
			gsl_vector_get(x, n) + gsl_ran_gaussian(r, step_size));
	}

	if(!gsl_finite(loglik_old))
		return 0;

	loglik_new = logdistr(data, x);
	if(loglik_new >= loglik_old)
		return 0;

	lik_ratio = exp(loglik_new - loglik_old);
	if(gsl_finite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio))
		return 0;

	gsl_vector_memcpy(x, old);
	return 0;
}
