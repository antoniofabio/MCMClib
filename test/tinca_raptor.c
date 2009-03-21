/**Test INCA RAPTOR algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt_q.h>
#include <inca_raptor.h>

#define N 100
#define DIM 1
#define K 2 /*number of regions*/
#define M 3 /*number of parallel chains*/
/*burn-in*/
#define T0 2
#define SF (2.38*2.38/(double) DIM)
#define CORRECTION 0.001

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_vector* x[M];
  for(int m=0; m<M; m++) {
    x[m] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(x[m], gsl_rng_uniform(rng));
  }

  gsl_matrix* Sigma_zero = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix* Sigma_hat[K];
  gsl_vector* mu_hat[K];
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], 0.0);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
  }
  gsl_vector* beta_hat = gsl_vector_alloc(K);
  gsl_vector_set_all(beta_hat, 1.0 / (double) K);

  mcmclib_inca* s = mcmclib_inca_raptor_alloc(rng,
					      dunif, NULL,
					      x, M, T0, Sigma_zero,
					      beta_hat, mu_hat, Sigma_hat);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++) {
    mcmclib_inca_update(s);
  }

  /*free memory*/
  for(int k=0; k<K; k++)
    gsl_matrix_free(Sigma_hat[k]);
  gsl_matrix_free(Sigma_zero);
  for(int m=0; m<M; m++)
    gsl_vector_free(x[m]);
  mcmclib_inca_raptor_free(s);

  return 0;
}
