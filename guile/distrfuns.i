%{
#include <mvnorm.h>
%}

/*Multivariate normal distribution*/
void mcmclib_mvnorm_iid(const gsl_rng* r, gsl_vector* out);
void mcmclib_mvnorm(const gsl_rng* r,
		    const gsl_matrix* sigma,
		    gsl_vector* out);
void mcmclib_mvnorm_chol(const gsl_rng* r,
	const gsl_matrix* sigma_chol,
	gsl_vector* out);
void mcmclib_mvnorm_precision(const gsl_rng* r,
			      const gsl_matrix* Psi,
			      gsl_vector* out);
void mcmclib_mvnorm_cholprec(const gsl_rng* r,
			     const gsl_matrix* Psi,
			     gsl_vector* out);

typedef struct {
} mcmclib_mvnorm_lpdf;
%extend mcmclib_mvnorm_lpdf {
  mcmclib_mvnorm_lpdf(gsl_vector* mean, double* vcov) {
    return mcmclib_mvnorm_lpdf_alloc(mean, vcov);
  }
  ~mcmclib_mvnorm_lpdf() {
    mcmclib_mvnorm_lpdf_free($self);
  }
}

double mcmclib_mvnorm_lpdf_compute(void* in_p, gsl_vector* x);
int mcmclib_mvnorm_lpdf_chol(mcmclib_mvnorm_lpdf* p);
double mcmclib_mvnorm_lpdf_compute_nochol(mcmclib_mvnorm_lpdf* p, gsl_vector* x);
void mcmclib_mvnorm_lpdf_inverse(mcmclib_mvnorm_lpdf* p);
double mcmclib_mvnorm_lpdf_compute_noinv(mcmclib_mvnorm_lpdf* p, gsl_vector* x);
double mcmclib_mvnorm_lpdf_noinv(gsl_vector* mu, gsl_matrix* iSigma, gsl_vector* x,
				 double ldet, gsl_vector* work1, gsl_vector* work2);
double mcmclib_mvnormzp_lpdf(const gsl_matrix* Psi, const gsl_vector* y);
