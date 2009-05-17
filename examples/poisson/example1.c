/**Adaptive Gaussian Random Walk on an MCAR model example*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>
#include <pois_model.h>

#define N 100000
#define THIN 10
#define T0 10000
#define V0 0.4

#define n 6
#define P 3

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
  gsl_matrix* X = gsl_matrix_alloc(n, P);
  gsl_matrix_set_zero(X);
  for(int i=0; i<P; i++)
    gsl_matrix_set(X, i, i, 1.0);
  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector_set_all(y, 3.0);
  mcmclib_pmodel_sampler* model = mcmclib_pmodel_sampler_alloc(X, y, NULL, rng, V0, T0);

  FILE* out = fopen("chain_beta.dat", "w");
  for(int i=0; i<N; i++) {
    if(((i+1) % THIN) == 0)
      gsl_vector_fprintf(out, mcmclib_pmodel_sampler_beta(model), "%f");
    mcmclib_pmodel_sampler_update(model);
  }
  fclose(out);

  mcmclib_pmodel_sampler_free(model);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_rng_free(rng);
}
