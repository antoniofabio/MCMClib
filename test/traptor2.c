/**Test RAPTOR: check custom alpha fun. setting*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <region_mixnorm.h>
#include <raptor.h>

static const double beta = 0.5;
static const double V[] = {1.0, 1.0};
static const double MU[] = {0.2, 0.8};

#define N 100
#define DIM 1
#define K 2
#define T0 50

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static double volume = 1.0;
static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= volume))
    return 0.0;
  return log(0.0);
}

static int check_alpha_ncalls = 0;
static double test_alpha_fun(void* ignore, mcmclib_raptor_gamma* g) {
  check_alpha_ncalls++;
  return 0.1;
}

static double test_alpha_fun_2(void* ignore, mcmclib_raptor_gamma* g) {
  return mcmclib_raptor_alpha_star_fun(g);
}

int main(int argc, char** argv) {
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);

  gsl_vector* mu_hat[K];
  gsl_matrix* Sigma_hat[K];
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], MU[k] * 0.5);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
  }
  gsl_vector* w_hat = gsl_vector_alloc(K);
  gsl_vector_set_all(w_hat, 1.0 / (double) K);

  mcmclib_amh* sampler = mcmclib_raptor_alloc(rng,
					      dunif, NULL,
					      x, T0, sigma_whole,
					      w_hat, mu_hat, Sigma_hat);

  /*check options setting*/
  mcmclib_raptor_set_alpha_fun(sampler, NULL, test_alpha_fun);
  for(int n=0; n<N; n++)
    mcmclib_amh_update(sampler);
  assert(check_alpha_ncalls == (N - T0));

  mcmclib_raptor_set_alpha_fun_identity(sampler);
  for(int n=0; n<N; n++)
    mcmclib_amh_update(sampler);

  mcmclib_raptor_set_alpha_fun(sampler, NULL, test_alpha_fun_2);
  for(int n=0; n<N; n++)
    mcmclib_amh_update(sampler);

  /*free memory*/
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_raptor_free(sampler);
  gsl_vector_free(w_hat);
  for(int k=0; k<K; k++) {
    gsl_vector_free(mu_hat[k]);
    gsl_matrix_free(Sigma_hat[k]);
  }
  gsl_rng_free(rng);

  return 0;
}
