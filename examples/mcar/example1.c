/**Adaptive Gaussian Random Walk on an MCAR model example*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>
#include <mcar_model.h>

#define N 1000
#define T0 50
#define V0 0.1

#define P 3
#define DIM 10

mcmclib_mcar_model* mcmclib_mcar_model_alloc(mcmclib_mcar_tilde_lpdf* m, gsl_vector* e);
void mcmclib_mcar_model_free(mcmclib_mcar_model* p);

double mcmclib_mcar_model_alpha1_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha1);
double mcmclib_mcar_model_alpha2_lpdf(mcmclib_mcar_model* p, gsl_vector* alpha2);
double mcmclib_mcar_model_sigma_lpdf(mcmclib_mcar_model* p, gsl_vector* sigma);
double mcmclib_mcar_model_Gamma_lpdf(mcmclib_mcar_model* p, gsl_vector* gamma);
double mcmclib_mcar_model_alphasigma_lpdf(mcmclib_mcar_model* p, gsl_vector* alphasigma);

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
  mcmclib_mcar_model* model = mcmclib_mcar_model_alloc(lpdf, e);

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  gsl_rng_free(r);
  mcmclib_mcar_model_free(model);
  gsl_vector_free(e);
  mcmclib_mcar_tilde_lpdf_free(lpdf);
}
