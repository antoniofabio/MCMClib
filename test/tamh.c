/**Test MH object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <amh.h>

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static double dtarget(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  if((x0 < -10.0))
    return log(0.5);
  return log(0.0);
}

static double qd(void* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

static void sampler(void* state, gsl_vector* x) {
  double* o = (double*) state;
  gsl_vector_set(x, 0, x0 + (*o));
}

static void update_gamma(void* in_p, gsl_vector* ignore) {
  int *p = (int*) in_p;
  (*p)++;
}

int main(int argc, char** argv) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, -0.5);
  double inc = 0.2;
  int gamma = 0;
  mcmclib_mh* mh = mcmclib_mh_alloc(r, dtarget, NULL,
				    x, qd, NULL,
				    sampler, &inc);
  mcmclib_amh* s = mcmclib_amh_alloc(mh, update_gamma, &gamma);

  for(int n=0; n<10; n++) {
    mcmclib_amh_update(s);
  }
  assert(gamma == 10);

  gsl_vector_free(x);
  gsl_rng_free(r);
  mcmclib_amh_free(s);
  mcmclib_mh_free(mh);

  return 0;
}
