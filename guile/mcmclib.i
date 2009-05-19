%module mcmclib

%include "typemaps.i"
%include "carrays.i"
%include "gsl.i"

%array_functions(double, doubleArray);
%array_functions(char *, stringArray);
%array_functions(int, intArray);

%{
#include <lpdf_iwishart.h>
#include <mcar_tilde.h>
#include <mcar_model.h>
#include <pois_model.h>
%}

mcmclib_iwishart_lpdf* mcmclib_iwishart_lpdf_alloc(gsl_matrix* Psi, int m);
void mcmclib_iwishart_lpdf_free(mcmclib_iwishart_lpdf* p);
double mcmclib_iwishart_lpdf_compute(void* p, gsl_vector* x);

double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x);
void mcmclib_mcar_tilde_lpdf_update_B_tilde(mcmclib_mcar_tilde_lpdf* p);
int mcmclib_mcar_tilde_lpdf_update_blocks(mcmclib_mcar_tilde_lpdf* p);
int mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p);
void mcmclib_mcar_tilde_lpdf_update_Gamma(mcmclib_mcar_tilde_lpdf* p);

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

  /*workspace memory*/
  gsl_matrix *Lambda_ij, *Gammai, *Block; /**< used in fun 'vcov_blockij' */
} mcmclib_mcar_tilde_lpdf;
%extend mcmclib_mcar_tilde_lpdf {
  mcmclib_mcar_tilde_lpdf(int p, gsl_matrix* M) {
    return mcmclib_mcar_tilde_lpdf_alloc(p, M);
  }
  ~mcmclib_mcar_tilde_lpdf() {
    mcmclib_mcar_tilde_lpdf_free($self);
  }
}

double mcmclib_mcar_model_alpha12sigma_lpdf(void* in_p, gsl_vector* alpha12sigma);
double mcmclib_mcar_model_alphasigma_lpdf(void* in_p, gsl_vector* alphasigma);
double mcmclib_mcar_model_phi_fcond(mcmclib_mcar_model* in_p, int i, gsl_vector* x);

typedef struct {
  mcmclib_mcar_tilde_lpdf* lpdf; /**< log-likelihood component */
  gsl_vector* e; /**< observed residuals*/
} mcmclib_mcar_model;
%extend mcmclib_mcar_model {
  mcmclib_mcar_model(mcmclib_mcar_tilde_lpdf* p, gsl_vector* e) {
    return mcmclib_mcar_model_alloc(p, e);
  }
  ~mcmclib_mcar_model() {
    mcmclib_mcar_model_free($self);
  }
}

int mcmclib_pois_model_set_prior_mean(mcmclib_pois_model* p, const gsl_vector* b0);
int mcmclib_pois_model_set_prior_var(mcmclib_pois_model* p, const gsl_matrix* B0);
int mcmclib_pois_model_set_offset(mcmclib_pois_model* p, const gsl_vector* offset);

double mcmclib_pois_model_llik(mcmclib_pois_model* p, gsl_vector* x);
double mcmclib_pois_model_lprior(mcmclib_pois_model* p, gsl_vector* x);
double mcmclib_pois_model_lpdf(void* in_p, gsl_vector* x);

typedef struct {
  gsl_vector* beta;

  /*internal*/
  gsl_matrix* X; /**< regression matrix */
  const gsl_vector* y; /**< observed values */
  const gsl_vector* offset; /**< (optional) offset on the mean */
  gsl_vector* b0; /**< prior beta mean */
  gsl_matrix* B0; /**< prior beta precision */
  gsl_vector* mu; /**< log-means (= X beta) */
  double ldet; /**< log-determinant of the inverse of B0 */
  gsl_vector *work1, *work2; /**< workspace memory */
} mcmclib_pois_model;

%extend mcmclib_pois_model {
  mcmclib_pois_model(gsl_matrix* X, gsl_vector* y) {
    return mcmclib_pois_model_alloc(X, y);
  }
  ~mcmclib_pois_model() {
    mcmclib_pois_model_free($self);
  }
}

typedef struct {
  mcmclib_pois_model* model;
  mcmclib_amh* sampler;
} mcmclib_pmodel_sampler;

mcmclib_pmodel_sampler* mcmclib_pmodel_sampler_alloc(const gsl_matrix* X,
						     const gsl_vector* y,
						     const gsl_vector* offset,
						     gsl_rng* rng,
						     double sigma0,
						     int burnin);

void mcmclib_pmodel_sampler_free(mcmclib_pmodel_sampler* p);
int mcmclib_pmodel_sampler_update(mcmclib_pmodel_sampler* p);
#define mcmclib_pmodel_sampler_beta(p) (p)->model->beta

%extend mcmclib_pmodel_sampler {
  mcmclib_pmodel_sampler(const gsl_matrix* X,
			 const gsl_vector* y,
			 const gsl_vector* offset,
			 gsl_rng* rng,
			 double sigma0,
			 int burnin) {
    return mcmclib_pmodel_sampler_alloc(X, y, offset, rng, sigma0, burnin);
  }
  ~mcmclib_pmodel_sampler() {
    mcmclib_pmodel_sampler_free($self);
  }
}
