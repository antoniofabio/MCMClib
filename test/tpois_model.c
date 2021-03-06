/**Test MCAR model */
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <pois_model.h>
#include "CuTest.h"

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define P 3
#define N 4

static gsl_vector* x;
static mcmclib_pois_model* mod;

static double llik(double sx) {
  gsl_vector_set_all(x, sx);
  return mcmclib_pois_model_llik(mod, x);
}

static double lprior(double sx) {
  gsl_vector_set_all(x, sx);
  return mcmclib_pois_model_lprior(mod, x);
}

static double lpdf(double sx) {
  gsl_vector_set_all(x, sx);
  return mcmclib_pois_model_lpdf(mod, x);
}

void Testpois_model(CuTest* tc) {
  gsl_matrix* X = gsl_matrix_alloc(N, P);
  gsl_matrix_set_zero(X);
  for(size_t i=0; i<P; i++)
    gsl_matrix_set(X, i, i, 1.0);
  gsl_vector* y = gsl_vector_alloc(N);
  gsl_vector_set_all(y, 2.0);
  mod = mcmclib_pois_model_alloc(X, y);

  CuAssertTrue(tc, mod->beta->size == P);
  x = gsl_vector_alloc(P);
  llik(0.0);
  llik(1.0);
  llik(2.0);

  lprior(1.0);
  CuAssertTrue(tc, check_dequal(lpdf(1.0), lprior(1.0) + llik(1.0)));

  gsl_vector* b0 = gsl_vector_alloc(P);
  gsl_vector_set_all(b0, 1.0);
  mcmclib_pois_model_set_prior_mean(mod, b0);
  gsl_vector_free(b0);
  gsl_matrix* B0 = gsl_matrix_alloc(P, P);
  gsl_matrix_set_identity(B0);
  mcmclib_pois_model_set_prior_var(mod, B0);
  gsl_matrix_free(B0);

  gsl_vector* offset = gsl_vector_alloc(N);
  gsl_vector_set_all(offset, -1.0);
  mcmclib_pois_model_set_offset(mod, offset);
  llik(1.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
  mcmclib_pmodel_sampler* ms = mcmclib_pmodel_sampler_alloc(X, y, offset, rng, 0.4, 100);
  mcmclib_pmodel_sampler_update(ms);

  mcmclib_pmodel_sampler_free(ms);
  gsl_rng_free(rng);
  mcmclib_pois_model_free(mod);
  gsl_vector_free(offset);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(x);
}
