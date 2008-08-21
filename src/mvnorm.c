#include <math.h>
#include "mvnorm.h"

#define PI 3.14159265

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

mvnorm_lpdf_p* mcmclib_mvnorm_lpdf_alloc(gsl_vector* mean, gsl_matrix* vcov) {
	int d = mean->size;
	mvnorm_lpdf_p* ans = (mvnorm_lpdf_p*) malloc(sizeof(mvnorm_lpdf_p));
	ans->mean = mean;
	ans->vcov = vcov;
	ans->rooti = gsl_matrix_alloc(d, d);
	ans->x_mu = gsl_vector_alloc(d);
	return ans;
}

void mcmclib_mvnorm_lpdf_free(mvnorm_lpdf_p* p) {
	gsl_matrix_free(p->rooti);
	gsl_vector_free(p->x_mu);
	free(p);
}

double mcmclib_mvnorm_lpdf(gsl_vector* x, void* in_p) {
	int d = x->size;
	mvnorm_lpdf_p* p = (mvnorm_lpdf_p*) in_p;
	gsl_matrix_memcpy(p->rooti, p->vcov);
	gsl_linalg_cholesky_decomp(p->rooti);
	gsl_vector* x_mu = p->x_mu;
	gsl_vector_memcpy(x_mu, x);
	gsl_vector_sub(x_mu, p->mean);
	/*read the following as: z = as.vector(t(rooti) %*% (x - mu))*/
	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, p->rooti, x_mu);

	/*read the following block as:
		-(length(x)/2) * log(2 * pi) - 0.5 * (z %*% z) + sum(log(diag(rooti))) */
	double ans = 0.0;
	gsl_blas_ddot(x_mu, x_mu, &ans);
	ans *= -0.5;
	ans -= log(2.0 * PI) * ((double) d) / 2.0;
	for(int i=0; i<d; i++)
		ans += log(gsl_matrix_get(p->rooti, i, i));

	return ans;
}
