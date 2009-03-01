/**Test MH object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mh.h>

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

int main(int argc, char** argv) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, -0.5);
  double inc = 0.2;
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, sampler, &inc, qd, NULL, &inc);
  mcmclib_mh* s = mcmclib_mh_alloc(r, dtarget, NULL, q, x);

  for(int n=0; n<10; n++) {
    mcmclib_mh_update(s);
  }
  assert(check_dequal(x0, 0.9));
  inc = -0.2;
  mcmclib_mh_update(s);
  assert(check_dequal(x0, 0.7));
  inc = -1.0;
  mcmclib_mh_update(s);
  assert(check_dequal(x0, 0.7));

  inc = -11.0;
  int count = 0;
  for(int n=0; n<10000; n++) {
    gsl_vector_set(x, 0, 0.0);
    mcmclib_mh_update(s);
    count += check_dequal(x0, inc);
  }
  /*jump should be accepted almost 50% of the times:*/
  assert(abs(count - 5000) < 100);

  gsl_vector_free(x);
  gsl_rng_free(r);
  mcmclib_mh_free(s);

  return 0;
}
