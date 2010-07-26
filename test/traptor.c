/**Test RAPTOR algorithm on dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <region_mixnorm.h>
#include "raptor.c"
#include "CuTest.h"

static const double beta = 0.5;
static const double V[] = {1.0, 1.0};
static const double MU[] = {0.2, 0.8};

#define N 100
#define DIM 1
#define K 2
#define T0 50
#define SF (2.38*2.38/(double) DIM)

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

#define TOL 1e-6

static double dunif(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

double fix(double in) {
  return (in + 0.001) * SF;
}

void Testraptor(CuTest* tc) {
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);

  gsl_vector* mu_hat[K], *mu[K];
  gsl_matrix* Sigma_hat[K], *Sigma[K];
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], MU[k] * 0.5);
    mu[k] = gsl_vector_alloc(DIM);
    gsl_vector_memcpy(mu[k], mu_hat[k]);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
    Sigma[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_memcpy(Sigma[k], Sigma_hat[k]);
  }
  gsl_vector* w_hat = gsl_vector_alloc(K);
  gsl_vector_set_all(w_hat, 1.0 / (double) K);
  gsl_vector* beta_v = gsl_vector_alloc(K);
  gsl_vector_memcpy(beta_v, w_hat);
  mcmclib_mixem_online* olem = mcmclib_mixem_online_alloc(mu, Sigma, beta_v, 0.6, T0);

  mcmclib_amh* sampler = mcmclib_raptor_alloc(rng,
					      dunif, NULL,
					      x, T0, sigma_whole,
					      w_hat, mu_hat, Sigma_hat);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++) {
    mcmclib_amh_update(sampler);
    mcmclib_mixem_online_update(olem, sampler->mh->x);
  }

  mcmclib_rapt_gamma* q = (mcmclib_rapt_gamma*) sampler->mh->q->gamma;
  raptor_gamma* g = (raptor_gamma*) q->which_region_data;
  CuAssertDblEquals(tc, v0(olem->beta), v0(g->beta_hat), TOL);
  CuAssertDblEquals(tc, v0(olem->mu[0]), v0(g->mu_hat[0]), TOL);
  CuAssertDblEquals(tc, m00(olem->Sigma[0]), m00(g->Sigma_hat[0]), TOL);

  /*check boundary function*/
  gsl_vector_set_all(x, -1.0);
  size_t rx = mcmclib_region_mixnorm_compute(x, g->pi_hat);
  CuAssertIntEquals(tc, 0, rx);
  gsl_vector_set_all(x, 1.0);
  rx = mcmclib_region_mixnorm_compute(x, g->pi_hat);
  CuAssertIntEquals(tc, 1, rx);

  /*check options setting*/
  mcmclib_raptor_set_alpha(sampler, 0.2);
  for(int n=0; n<N; n++)
    mcmclib_amh_update(sampler);
  mcmclib_raptor_set_sf(sampler, 0.1);
  for(int n=0; n<N; n++)
    mcmclib_amh_update(sampler);

  /*free memory*/
  mcmclib_mixem_online_free(olem);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_amh_free(sampler);
  gsl_rng_free(rng);
  for(int k=0; k<K; k++) {
    gsl_vector_free(mu_hat[k]);
    gsl_vector_free(mu[k]);
    gsl_matrix_free(Sigma_hat[k]);
    gsl_matrix_free(Sigma[k]);
  }
  gsl_vector_free(w_hat);
  gsl_vector_free(beta_v);
}
