#include <gsl/gsl_math.h>
#include "metropolis.h"
#include "gauss_mrw.h"
#include "mvnorm.h"

mcmclib_gauss_mrw* mcmclib_gauss_mrw_alloc(gsl_rng* r,
	distrfun_p logdistr, void* logdistr_data, gsl_vector* start_x, const gsl_matrix* sigma_prop) {
	int dim = start_x->size;
	mcmclib_gauss_mrw* ans = (mcmclib_gauss_mrw*) malloc(sizeof(mcmclib_gauss_mrw));
	ans->r = r;
	ans->logdistr = logdistr;
	ans->logdistr_data = logdistr_data;
	ans->current_x = start_x;
	ans->old = gsl_vector_alloc(dim);
	ans->sigma_prop = gsl_matrix_alloc(dim, dim);
	gsl_matrix_memcpy(ans->sigma_prop, sigma_prop);
	return ans;
}

void mcmclib_gauss_mrw_free(mcmclib_gauss_mrw* p) {
	gsl_matrix_free(p->sigma_prop);
	gsl_vector_free(p->old);
	free(p);
}

int mcmclib_gauss_mrw_update(mcmclib_gauss_mrw* p) {
	gsl_vector_memcpy(p->old, p->current_x);
	mcmclib_mvnorm(p->r, p->sigma_prop, p->current_x);
	gsl_vector_add(p->current_x, p->old);
	int ans = mcmclib_metropolis_symmetric_step(p->r,
		p->old, p->current_x, p->logdistr, p->logdistr_data);
	return ans;
}
