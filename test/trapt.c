/**Test base RAPT algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mvnorm.h>
#include <rapt.h>

#define N 1000
#define DIM 1
#define K 2
/*burn-in*/
#define T0 100

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static int which_region(gsl_vector* x, void* ignore) {
  return gsl_vector_get(x, 0) < 0 ? 0 : 1;
}

int main(int argc, char** argv) {
  /*setup target distrib.*/
  gsl_vector* mu = gsl_vector_alloc(DIM);
  gsl_vector_set_all(mu, 0.0);
  gsl_matrix* Sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma);
  mcmclib_mvnorm_lpdf* pi = mcmclib_mvnorm_lpdf_alloc(mu, Sigma->data);

  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);
  gsl_matrix* sigma_local[K];
  for(int k=0; k<K; k++) {
    sigma_local[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(sigma_local[k]);
  }

  mcmclib_rapt* sampler = mcmclib_rapt_alloc(rng,
					     mcmclib_mvnorm_lpdf_compute, pi, /*target distrib.*/
					     x, T0,
					     sigma_whole, K, sigma_local,
					     which_region, NULL);

  /*Main MCMC loop*/
  double sum_x = 0.0;
  double sum_x2 = 0.0;
  for(int n=0; n<N; n++) {
    mcmclib_rapt_update(sampler);
    sum_x += gsl_vector_get(x, 0);
    mcmclib_rapt_update_proposals(sampler);
    sum_x2 += pow(gsl_vector_get(x, 0), 2);
  }
  /*check results*/
  assert(check_dequal(sum_x, -57.336886));
  assert(check_dequal(sum_x2, 912.202062));

  /*free memory*/
  for(int k=0; k<K; k++)
    gsl_matrix_free(sigma_local[k]);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_rapt_free(sampler);
  mcmclib_mvnorm_lpdf_free(pi);
  gsl_vector_free(mu);
  gsl_matrix_free(Sigma);

  return 0;
}
