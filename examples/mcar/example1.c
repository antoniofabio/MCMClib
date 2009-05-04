/**Adaptive Gaussian Random Walk on an MCAR model example*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>
#include <mcar_tilde.h>
#include <mcar_model.h>

#define N 10
#define T0 50
#define V0 0.4

#define P 3
#define DIM 2

gsl_vector *alpha12sigma, *alphasigmag;
gsl_rng* rng;
mcmclib_mcar_tilde_lpdf* lpdf;
mcmclib_mcar_model* model;
mcmclib_amh* sampler[2];

void init_chains() {
  rng = gsl_rng_alloc(gsl_rng_default);

  alpha12sigma = lpdf->alpha12sigma;
  gsl_vector_set_all(alpha12sigma, -1.0);
  gsl_matrix* Sigma0 = gsl_matrix_alloc(P*P, P*P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / ((double)(P * P)));
  sampler[0] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alpha12sigma_lpdf,
				      model, alpha12sigma, Sigma0, T0);
  gsl_matrix_free(Sigma0);

  alphasigmag = lpdf->alphasigmag;
  gsl_vector_set_all(alphasigmag, -1.0);
  Sigma0 = gsl_matrix_alloc(P*(P-1)/2 + P, P*(P-1)/2 + P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / (double) (P*(P-1)/2 + P));
  sampler[1] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alphasigma_lpdf,
				      model, alphasigmag, Sigma0, T0);
  gsl_matrix_free(Sigma0);
}

void free_chains() {
  for(int i=0; i<2; i++)
    mcmclib_gauss_am_free(sampler[i]);
  gsl_rng_free(rng);
}

int main(int argc, char** argv) {
  /* model setup */
  gsl_matrix* W = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_zero(W);
  for(int i=0; i<(DIM-1); i++) {
    gsl_matrix_set(W, i, i+1, 1.0);
    gsl_matrix_set(W, i+1, i, 1.0);
  }
  lpdf = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  gsl_matrix_free(W);
  gsl_vector* e = gsl_vector_alloc(P * DIM);
  gsl_vector_set_zero(e);
  model = mcmclib_mcar_model_alloc(lpdf, e);

  init_chains();
  for(int i=0; i<N; i++) {
    for(int j=0; j<2; j++) {
      mcmclib_amh_update(sampler[j]);
    }
  }

  free_chains();
  mcmclib_mcar_model_free(model);
  gsl_vector_free(e);
  mcmclib_mcar_tilde_lpdf_free(lpdf);
}