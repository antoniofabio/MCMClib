/**Adaptive Gaussian Random Walk and RAPTOR on a Poisson model
   with MCAR effects example*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <mvnorm.h>
#include <mcar_tilde.h>
#include <mcar_model.h>
#include <pois_model.h>
#include <gauss_am.h>
#include <raptor.h>

/*P=3, DIM=95: 0.21009 secs per iteration
  P=6, DIM=10: 0.00396 secs per iteration
  P=3, DIM=10: 0.00085 secs per iteration
  P=3, DIM=5:  0.00016 secs per iteration
*/
#define N 50000
#define THIN 10
#define T0 1000
#define V0 0.4

#define P 3
#define DIM 95

gsl_rng* rng;
gsl_vector *alpha12sigma, *alphasigmag;
mcmclib_mcar_tilde_lpdf* mcar_lpdf;
mcmclib_mcar_model* mcar_model;
gsl_vector* mcar_phi;
gsl_vector* y;
gsl_matrix* X;
mcmclib_pmodel_sampler* model;
mcmclib_amh* sampler[3];

void init_chains() {
  alpha12sigma = mcar_lpdf->alpha12sigma;
  gsl_vector_set_all(alpha12sigma, -1.0);
  gsl_matrix* Sigma0 = gsl_matrix_alloc(P*P, P*P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / ((double)(P * P)));
  sampler[0] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alpha12sigma_lpdf,
				      mcar_model, alpha12sigma, Sigma0, T0);
  gsl_matrix_free(Sigma0);

  alphasigmag = mcar_lpdf->alphasigmag;
  gsl_vector_set_all(alphasigmag, -1.0);
  Sigma0 = gsl_matrix_alloc(P*(P-1)/2 + P, P*(P-1)/2 + P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / (double) (P*(P-1)/2 + P));
  sampler[1] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alphasigma_lpdf,
				      mcar_model, alphasigmag, Sigma0, T0);
  gsl_matrix_free(Sigma0);

  model = mcmclib_pmodel_sampler_alloc(X, y, mcar_phi, rng, 1e-3, T0);
  gsl_vector_set_all(model->model->beta, 0.0);
  gsl_vector_set(model->model->beta, 0, -1.0);
  gsl_vector_set(model->model->beta, 2, 1.0);
  sampler[2] = model->sampler;
}

void free_chains() {
  for(int i=0; i<2; i++)
    mcmclib_gauss_am_free(sampler[i]);
  mcmclib_pmodel_sampler_free(model);
  gsl_rng_free(rng);
}

int main(int argc, char** argv) {
  rng = gsl_rng_alloc(gsl_rng_default);

  /* model setup */
  gsl_matrix* W = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_zero(W);
  for(int i=0; i<(DIM-1); i++) {
    gsl_matrix_set(W, i, i+1, 1.0);
    gsl_matrix_set(W, i+1, i, 1.0);
  }
  mcar_lpdf = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  gsl_matrix_free(W);
  mcar_phi = gsl_vector_alloc(P * DIM);
  gsl_vector_set_zero(mcar_phi);
  mcar_model = mcmclib_mcar_model_alloc(mcar_lpdf, mcar_phi);

  y = gsl_vector_alloc(P * DIM);
  FILE* in_y = fopen("y.dat", "r");
  gsl_vector_fscanf(in_y, y);
  fclose(in_y);
  X = gsl_matrix_alloc(P * DIM, P);
  gsl_matrix_set_zero(X);
  for(int d=0; d<DIM; d++)
    for(int i=0; i<P; i++)
      gsl_matrix_set(X, d*P + i, i, 1.0);

  init_chains();
  FILE* out_a12s = fopen("chain_alpha12sigma.dat", "w");
  FILE* out_as = fopen("chain_alphasigma.dat", "w");
  FILE* out_beta = fopen("chain_beta.dat", "w");
  FILE* out_phi = fopen("chain_phi.dat", "w");
  FILE* out_lpdf = fopen("chain_lpdf.dat", "w");
  for(int i=0; i<N; i++) {
    if (( (i+1) % THIN ) == 0) {
      gsl_vector_fprintf(out_a12s, alpha12sigma, "%f");
      gsl_vector_fprintf(out_as, alphasigmag, "%f");
      gsl_vector_fprintf(out_beta, sampler[2]->mh->x, "%f");
      fprintf(out_lpdf, "%f\n", mcmclib_pois_model_lpdf(model->model, sampler[2]->mh->x));
      gsl_vector_fprintf(out_phi, mcar_phi, "%f");
      fflush(out_beta);
      fflush(out_lpdf);
      fflush(out_phi);
    }
    for(int j=0; j<2; j++) {
      mcmclib_amh_update(sampler[j]);
      assert(gsl_finite(sampler[j]->mh->logdistr_old));
    }
    mcmclib_mcar_tilde_lpdf_update_vcov(mcar_lpdf);
    gsl_matrix* tmp = gsl_matrix_alloc(P*DIM, P*DIM);
    gsl_matrix_memcpy(tmp, mcar_lpdf->vcov);
    gsl_linalg_cholesky_decomp(tmp);
    mcmclib_mvnorm_precision(rng, tmp, mcar_phi);
    gsl_matrix_free(tmp);
    mcmclib_pmodel_sampler_update(model);
  }
  fclose(out_a12s);
  fclose(out_as);
  fclose(out_beta);
  fclose(out_phi);
  fclose(out_lpdf);

  free_chains();
  mcmclib_mcar_model_free(mcar_model);
  gsl_vector_free(mcar_phi);
  mcmclib_mcar_tilde_lpdf_free(mcar_lpdf);
}
