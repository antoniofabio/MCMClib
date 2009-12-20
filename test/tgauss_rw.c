/**Test GAUSS-RW algorithm on std. norm. target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_rw.h>

#define N 10000
#define DIM 1

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

static double dtarget(void* ignore, gsl_vector* x) {
  return log(gsl_ran_gaussian_pdf(gsl_vector_get(x, 0), 1.0));
}

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, gsl_rng_uniform(rng));

  double mean = 0.0;
  double variance = 0.0;
  double kurt = 0.0;

  mcmclib_mh* s = mcmclib_gauss_rw_alloc(rng,
					 dtarget, NULL, /*target distrib.*/
					 x, 1.0);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++) {
    mcmclib_mh_update(s);

    mean += x0;
    variance += x0 * x0;
    kurt += x0 * x0 * x0 * x0;
  }

  /*compute mean and variance*/
  mean /= (double) N;
  variance = variance / ((double) N) - (mean * mean);
  kurt /= (double) N;

  //assert(check_dequal(0.066613, mean));
  //assert(check_dequal(0.999890, variance));
  //assert(check_dequal(2.985794, kurt));

  /*free memory*/
  gsl_vector_free(x);
  mcmclib_amh_free(s);
  gsl_rng_free(rng);

  return 0;
}
