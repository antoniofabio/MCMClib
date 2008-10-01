#include <math.h>
#include <gsl/gsl_math.h>
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

mcmclib_mvnorm_lpdf* mcmclib_mvnorm_lpdf_alloc(gsl_vector* mean, double* vcov) {
	int d = mean->size;
	mcmclib_mvnorm_lpdf* ans = (mcmclib_mvnorm_lpdf*) malloc(sizeof(mcmclib_mvnorm_lpdf));
	ans->mean = mean;
	ans->vcov = vcov;
	ans->rooti = gsl_matrix_alloc(d, d);
	ans->x_mu = gsl_vector_alloc(d);
	ans->mahal = gsl_vector_alloc(d);
	return ans;
}

void mcmclib_mvnorm_lpdf_free(mcmclib_mvnorm_lpdf* p) {
	gsl_matrix_free(p->rooti);
	gsl_vector_free(p->x_mu);
	gsl_vector_free(p->mahal);
	free(p);
}

double mcmclib_mvnorm_lpdf_compute(void* in_p, gsl_vector* x) {
	int d = x->size;
	mcmclib_mvnorm_lpdf* p = (mcmclib_mvnorm_lpdf*) in_p;
	gsl_matrix_view mv = gsl_matrix_view_array(p->vcov, d, d);
	gsl_matrix* vcov = &(mv.matrix);

	/*compute cholesky decomposition of var/cov matrix*/
	gsl_matrix_memcpy(p->rooti, vcov);
	gsl_linalg_cholesky_decomp(p->rooti);

	/*compute mahlanobis distance between 'x' and 'mu'*/
	gsl_vector* x_mu = p->x_mu;
	gsl_vector_memcpy(x_mu, x);
	gsl_vector_sub(x_mu, p->mean);
	gsl_vector_memcpy(p->mahal, x_mu);
	gsl_linalg_cholesky_svx(p->rooti, x_mu);

	/*compute log-density as:
		-0.5 * (mahaldist + log(2*pi)*d + logdet) */
	double ans = 0.0;
	gsl_blas_ddot(p->mahal, x_mu, &ans);
	ans += log(2.0 * M_PI) * ((double) d);
	for(int i=0; i<d; i++)
		ans += log(gsl_matrix_get(p->rooti, i, i)) * 2.0;
	ans *= -0.5;

	return ans;
}
