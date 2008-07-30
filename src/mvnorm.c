#include "mvnorm.h"

void mcmclib_mvnorm(const gsl_rng* r,
	const gsl_matrix* sigma,
	gsl_vector* out) {
	int d = sigma->size1;

	gsl_matrix* sigma_chol = gsl_matrix_alloc(d, d);
	gsl_matrix_memcpy(sigma_chol, sigma);
	gsl_linalg_cholesky_decomp(sigma_chol);

	mcmclib_mvnorm_chol(r, sigma_chol, out);
	
	gsl_matrix_free(sigma_chol);
}

void mcmclib_mvnorm_chol(const gsl_rng* r,
	const gsl_matrix* sigma_chol,
	gsl_vector* out) {
	/*generate d iid values*/
	for(int n= (sigma_chol->size1 - 1); n>=0; n--)
		gsl_vector_set(out, n, gsl_ran_gaussian(r, 1.0));

	/*rotate them according to the cholesky 'square root' of sigma*/
	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, sigma_chol, out);
}
