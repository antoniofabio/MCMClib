/**Test GAUSS-AM algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <scam.h>

#define N 5
#define DIM 2
#define T0 2

#define v0(x) gsl_vector_get(x, 0)
#define v1(x) gsl_vector_get(x, 1)
#define x0 v0(x)
#define x1 v1(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static int counter1 = 0;

static double dunif_j(void* ignore, int j, double xj) {
	counter1 ++;
	double min = j ? 1.0 : 0.0;
	double max = j ? 2.0 : 1.0;
  if((xj >= min) && (xj <= max)) {
    return log(1.0);
	}
  return log(0.0);
}

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_vector* x = gsl_vector_alloc(DIM);
	gsl_vector_set(x, 0, gsl_rng_uniform(rng));
	gsl_vector_set(x, 1, 1.0 + gsl_rng_uniform(rng));

  mcmclib_scam* s = mcmclib_scam_alloc(rng, dunif_j, NULL,
																			 x, 0.1, T0);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++)
    mcmclib_scam_update(s);

  /*check results*/
	assert(counter1 == (N * 2 * DIM)); /*check that logdistr. has been called the
																			 right number of times*/
	assert((gsl_vector_get(x, 0) >= 0.0) && (gsl_vector_get(x, 0) <= 1.0));
	assert((gsl_vector_get(x, 1) >= 1.0) && (gsl_vector_get(x, 1) <= 2.0));

  /*free memory*/
  mcmclib_scam_free(s);
  gsl_vector_free(x);
	gsl_rng_free(rng);

  return 0;
}
