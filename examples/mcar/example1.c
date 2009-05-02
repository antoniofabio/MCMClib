/**Adaptive Gaussian Random Walk on an MCAR model example*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>
#include <mcar_model.h>

#define N 1000
#define T0 50
#define V0 4.0

#define P 3
#define DIM 5
#define ALPHAP P*(P-1)/2

gsl_vector *alpha1, *alpha2, *sigma, *Gammav;
gsl_rng* rng;
mcmclib_mcar_model* model;
mcmclib_amh* sampler[4];

void init_chains() {
  rng = gsl_rng_alloc(gsl_rng_default);

  alpha1 = gsl_vector_alloc(ALPHAP);
  gsl_vector_set_zero(alpha1);
  gsl_matrix* Sigma0 = gsl_matrix_alloc(ALPHAP, ALPHAP);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / ALPHAP);
  sampler[0] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alpha1_lpdf,
				      model, alpha1, Sigma0, T0);

  alpha2 = gsl_vector_alloc(ALPHAP);
  gsl_vector_set_zero(alpha2);
  sampler[1] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alpha2_lpdf,
				      model, alpha2, Sigma0, T0);
  gsl_matrix_free(Sigma0);

  sigma = gsl_vector_alloc(P);
  gsl_vector_set_zero(sigma);
  Sigma0 = gsl_matrix_alloc(P, P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / P);
  sampler[2] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_sigma_lpdf,
				      model, sigma, Sigma0, T0);
  gsl_matrix_free(Sigma0);

  Gammav = gsl_vector_alloc(ALPHAP + P);
  gsl_vector_set_zero(Gammav);
  Sigma0 = gsl_matrix_alloc(ALPHAP + P, ALPHAP + P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / (ALPHAP + P));
  sampler[3] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alphasigma_lpdf,
				      model, Gammav, Sigma0, T0);
  gsl_matrix_free(Sigma0);
}

void free_chains() {
  for(int i=0; i<4; i++)
    mcmclib_gauss_am_free(sampler[i]);
  gsl_rng_free(rng);
  gsl_vector_free(alpha1);
  gsl_vector_free(alpha2);
  gsl_vector_free(sigma);
  gsl_vector_free(Gammav);
}

int main(int argc, char** argv) {
  /* model setup */
  gsl_matrix* W = gsl_matrix_alloc(DIM, DIM);
  FILE* W_file = fopen("W.txt", "r");
  if(!W_file) {
    printf("file 'W.txt' not found\n");
    exit(1);
  } 
  gsl_matrix_fscanf(W_file, W);
  fclose(W_file);
  mcmclib_mcar_tilde_lpdf* lpdf = mcmclib_mcar_tilde_lpdf_alloc(P, DIM, W);
  gsl_matrix_free(W);
  gsl_vector* e = gsl_vector_alloc(P * DIM);
  model = mcmclib_mcar_model_alloc(lpdf, e);

  init_chains();
  for(int i=0; i<N; i++) {
    for(int j=0; j<4; j++) {
      mcmclib_amh_update(sampler[j]);
    }
  }

  free_chains();
  mcmclib_mcar_model_free(model);
  gsl_vector_free(e);
  mcmclib_mcar_tilde_lpdf_free(lpdf);
}
