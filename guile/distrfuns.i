%{
#include <mvnorm.h>
#include <mcar_tilde.h>
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

/*MCAR distribution*/
typedef struct {
  int p; /**< dimension */
  int n; /**< number of points */

  gsl_matrix* B_tilde; /**< variance par. matrix (p x p) */
  gsl_vector* alpha12sigma; /**< Givens angles and sing. values repr. of
			       B_tilde */
  gsl_matrix* Gamma; /**< 'variance of variance' par. matrix (p x p) */
  gsl_vector *alphasigmag; /**< Givens angles and eigenv. repr. of Gamma */

  gsl_matrix* M; /**< adiancency matrix (n x n)*/
  gsl_vector* m; /**< adiancency weights (n)*/

  gsl_matrix* vcov; /**< precision matrix */
} mcmclib_mcar_tilde_lpdf;

%extend mcmclib_mcar_tilde_lpdf {
  mcmclib_mcar_tilde_lpdf(int p, gsl_matrix* M) {
    return mcmclib_mcar_tilde_lpdf_alloc(p, M);
  }
  ~mcmclib_mcar_tilde_lpdf() {
    mcmclib_mcar_tilde_lpdf_free($self);
  }
}
%callback("%s_cb");
double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x);
%nocallback;

void mcmclib_mcar_tilde_lpdf_update_B_tilde(mcmclib_mcar_tilde_lpdf* p);
int mcmclib_mcar_tilde_lpdf_update_blocks(mcmclib_mcar_tilde_lpdf* p);
int mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p);
void mcmclib_mcar_tilde_lpdf_update_Gamma(mcmclib_mcar_tilde_lpdf* p);
