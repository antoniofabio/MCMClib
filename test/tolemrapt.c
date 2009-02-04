/**Test OLEM-RAPT algorithm on dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <region_mixnorm.h>
#include <olem_rapt.h>

static const double beta = 0.5;
static const double V[] = {1.0, 1.0};
static const double MU[] = {0.2, 0.8};

#define N 1000
#define DIM 1
#define K 2
#define T0 100
#define SF (2.38*2.38/(double) DIM)

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

#define TOL 1e-6
int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

double fix(double in, double correction) {
  return (in + correction) * SF;
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

  mcmclib_olemrapt* sampler = mcmclib_olemrapt_alloc(rng,
						     dunif, NULL,
						     x, T0, sigma_whole,
						     w_hat, mu_hat, Sigma_hat);

  /*Main MCMC loop*/
  double mean = 0.0;
  double var = 0.0;
  for(int n=0; n<N; n++) {
    mcmclib_olemrapt_update(sampler);
    mcmclib_olemrapt_update_proposals(sampler);
    mean += v0(x);
    var += pow(v0(x), 2);
  }
  mean /= (double) N;
  var /= (double) N;

  /*check boundary function*/
  gsl_vector_set_all(x, -1.0);
  int rx = mcmclib_region_mixnorm_compute(x, sampler->pi_hat);
  assert(rx == 0);
  gsl_vector_set_all(x, 1.0);
  rx = mcmclib_region_mixnorm_compute(x, sampler->pi_hat);
  assert(rx == 1);

  /*free memory*/
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_olemrapt_free(sampler);

  return 0;
}
