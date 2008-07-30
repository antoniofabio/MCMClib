#include "common.h"

int mcmclib_gauss_am(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	gsl_matrix* sigma_zero, int t0,
	gsl_matrix* cov, gsl_vector* mean, int* t) {
	//TODO
	return 0;
}
