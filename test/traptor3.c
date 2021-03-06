/**Test RAPTOR algorithm ergodicity on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <raptor.h>
#include "CuTest.h"

static double beta = 0.8;
static const double V[] = {1.0, 1.0};
static const double MU[] = {-0.5, 0.5};

#define N 1000
#define DIM 1
#define K 2
#define T0 500
#define SF 0.5

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static double dber(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  if(fabs(x0) > 2.0)
    return log(0.0);
  else if(x0 >= 0.0)
    return log(beta);
  else
    return log(1.0-beta);
}

void Testraptor3(CuTest* tc) {
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);

  gsl_vector* mu_hat[K];
  gsl_matrix* Sigma_hat[K];
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], MU[k]);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
  }
  gsl_vector* w_hat = gsl_vector_alloc(K);
  gsl_vector_set_all(w_hat, 1.0 / (double) K);
  mcmclib_amh* sampler = mcmclib_raptor_alloc(rng,
					      dber, NULL,
					      x, T0, sigma_whole,
					      w_hat, mu_hat, Sigma_hat);
  mcmclib_raptor_set_alpha(sampler, 0.3);

  /*Main MCMC loop*/
  for(int n=0; n<T0; n++) {
    mcmclib_amh_update(sampler);
  }
  int n1 = 0;
  int nacc = 0.0;
  for(int n=0; n<N; n++) {
    nacc += mcmclib_amh_update(sampler);
    n1 += (x0 >= 0.0) ? 1 : 0;
  }

  CuAssertIntEquals(tc, 794, n1);
  CuAssertIntEquals(tc, 458, nacc);

  /*free memory*/
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_amh_free(sampler);
  gsl_rng_free(rng);
  for(int k=0; k<K; k++) {
    gsl_vector_free(mu_hat[k]);
    gsl_matrix_free(Sigma_hat[k]);
  }
  gsl_vector_free(w_hat);
}
