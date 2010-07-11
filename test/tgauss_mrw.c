/**Test GAUSS-MRW algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_mrw.h>

#define N 1000
#define DIM 1

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static double dtarget(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  return log(gsl_ran_gaussian_pdf(x0, 1.0));
}

int main() {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, gsl_rng_uniform(rng));

  gsl_matrix* sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma);

  double mean = 0.0;
  double variance = 0.0;

  mcmclib_mh* s = mcmclib_gauss_mrw_alloc(rng,
					  dtarget, NULL, /*target distrib.*/
					  x, sigma);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++) {
    mcmclib_mh_update(s);
    mean += x0;
    variance += x0 * x0;
  }

  /*compute mean and variance*/
  mean /= (double) N;
  variance = variance / ((double) N) - (mean * mean);

  assert(check_dequal(mean, 0.072261));
  assert(check_dequal(variance, 0.927399));

  /*free memory*/
  gsl_matrix_free(sigma);
  gsl_vector_free(x);
  mcmclib_mh_free(s);
  gsl_rng_free(rng);

  return 0;
}
