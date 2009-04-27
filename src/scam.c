#include <assert.h>
#include "scam.h"
#include "gauss_am.h"

mcmclib_scam* mcmclib_scam_alloc(gsl_rng* r,
																 distrfun_i_p logdistr, void* logdistr_data,
																 gsl_vector* x,
																 double sigma_zero, int t0) {
	assert(sigma_zero > 0);
	mcmclib_scam* a = (mcmclib_scam*) malloc(sizeof(mcmclib_scam));
	a->x_full = x;
	int d = x->size;
	gsl_matrix* S0 = gsl_matrix_alloc(1, 1);
	gsl_matrix_set(S0, 0, 0, sigma_zero);
	a->logdistr = logdistr;
	a->xi = gsl_vector_alloc(1);
	a->x_smp = (mcmclib_amh**) malloc(d * sizeof(mcmclib_amh*));
	for(int i=0; i<d; i++)
		a->x_smp[i] = mcmclib_gauss_am_alloc(r, mcmclib_scam_logdistr, a,
																				 a->xi, S0, t0);
	gsl_matrix_free(S0);
	return a;
}

void mcmclib_scam_free(mcmclib_scam* p) {
	for(int i=0; i < p->x_full->size; i++)
		mcmclib_gauss_am_free(p->x_smp[i]);
	free(p->x_smp);
	gsl_vector_free(p->xi);
	free(p);
}

double mcmclib_scam_logdistr(void* data, gsl_vector* x) {
	mcmclib_scam* p = (mcmclib_scam*) data;
	return p->logdistr(p->logdistr_data, p->curr_index,
										 gsl_vector_get(x, 0));
}

void mcmclib_scam_update(mcmclib_scam* p) {
	for(int i=0; i< p->x_full->size; i++) {
		p->curr_index = i;
		gsl_vector_set(p->xi, 0, gsl_vector_get(p->x_full, i));
		mcmclib_amh_update(p->x_smp[i]);
		gsl_vector_set(p->x_full, i, gsl_vector_get(p->xi, 0));
	}
}
