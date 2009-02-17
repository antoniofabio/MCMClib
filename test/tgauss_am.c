/**Test GAUSS-AM algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>

#define N 1000
#define DIM 1
/*burn-in*/
#define T0 100
#define SF (2.38*2.38/(double) DIM)

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

static double fix(double in, double correction) {
  return (in + correction) * SF;
}

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, gsl_rng_uniform(rng));

  gsl_matrix* sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma);

  double mean = 0.0;
  double variance = 0.0;

  mcmclib_gauss_am* s = mcmclib_gauss_am_alloc(rng,
					       dunif, NULL, /*target distrib.*/
					       x, sigma, T0);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++) {
    mcmclib_gauss_am_update(s);

    mean += x0;
    variance += x0 * x0;
  }

  /*compute mean and variance*/
  mean /= (double) N;
  variance = variance / ((double) N) - (mean * mean);

  /*check results*/
  assert(s->amh->n == N);

  assert(check_dequal(mean, v0(s->mean)));
  assert(check_dequal(variance, m00(s->cov)));
  double eps = m00(s->Sigma_eps);
  mcmclib_gauss_mrw* mrw = (mcmclib_gauss_mrw*) s->mrw;
  assert(check_dequal(fix(variance, eps), m00(mrw->sigma_prop)));

  /*free memory*/
  gsl_matrix_free(sigma);
  gsl_vector_free(x);
  mcmclib_gauss_am_free(s);

  return 0;
}
