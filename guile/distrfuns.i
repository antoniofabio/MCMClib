%{
#include <mvnorm.h>
#include <mixnorm.h>
#include <mcar_tilde.h>
#include <mcar_model.h>
#include <pois_model.h>
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

%array_functions(mcmclib_mvnorm_lpdf*, mvnormArray);

typedef struct {
} mcmclib_mixnorm_lpdf;
%extend mcmclib_mixnorm_lpdf {
  mcmclib_mixnorm_lpdf(gsl_vector* w, mcmclib_mvnorm_lpdf** pis) {
    return mcmclib_mixnorm_lpdf_alloc(w, pis);
  }
  ~mcmclib_mixnorm_lpdf() {
    mcmclib_mixnorm_lpdf_free($self);
  }
}
%callback("%s_cb");
double mcmclib_mixnorm_lpdf_compute(void* p, gsl_vector* x);
%nocallback;

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

typedef struct {
  mcmclib_mcar_tilde_lpdf* lpdf; /**< log-likelihood component */
  gsl_vector* e; /**< observed residuals*/
} mcmclib_mcar_model;

%extend mcmclib_mcar_model {
  mcmclib_mcar_model(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e) {
    return mcmclib_mcar_model_alloc(m, e);
  }
  ~mcmclib_mcar_model() {
    mcmclib_mcar_model_free($self);
  }
}

%callback("%s_cb");
double mcmclib_mcar_model_alpha12sigma_lpdf(void* in_p, gsl_vector* alpha12sigma);
double mcmclib_mcar_model_alphasigma_lpdf(void* in_p, gsl_vector* alphasigma);
%nocallback;
double mcmclib_mcar_model_phi_fcond(mcmclib_mcar_model* in_p, int i, gsl_vector* x);

/*Poisson model sampler*/
typedef struct {
  mcmclib_pois_model* model;
  mcmclib_amh* sampler;
} mcmclib_pmodel_sampler;

%extend mcmclib_pmodel_sampler {
  mcmclib_pmodel_sampler(const gsl_matrix* X,
			 const gsl_vector* y,
			 const gsl_vector* offset,
			 gsl_rng* rng,
			 double sigma0,
			 int burnin) {
    return mcmclib_pmodel_sampler_alloc(X,
					y,
					offset,
					rng,
					sigma0,
					burnin);
  }
  ~mcmclib_pmodel_sampler() {
    mcmclib_pmodel_sampler_free($self);
  }
}
int mcmclib_pmodel_sampler_update(mcmclib_pmodel_sampler* p);
gsl_vector* mcmclib_pmodel_sampler_beta(mcmclib_pmodel_sampler* p);
