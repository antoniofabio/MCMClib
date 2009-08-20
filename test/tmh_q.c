/**Test MH-Q object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mh_q.h>

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

static double qd(void* in_inc, gsl_vector* x, gsl_vector* y) {
  double* inc = (double*) in_inc;
  return (*inc) * (v0(x) - v0(y));
}

static void sampler(mcmclib_mh_q* q, gsl_vector* x) {
  double* o = (double*) q->gamma;
  gsl_vector_set(x, 0, x0 + (*o));
}

int main(int argc, char** argv) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, 2.0);
  gsl_vector* y = gsl_vector_alloc(1);
  gsl_vector_set(y, 0, 4.0);
  double inc = 1.0;
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, sampler, &inc, qd, &inc, &inc);

  assert(mcmclib_mh_q_logd(q, x, y) == qd(&inc, x, y));
  assert(mcmclib_mh_q_ratio_offset(q, x, y) == qd(&inc, y, x) - qd(&inc, x, y));
  mcmclib_mh_q_sample(q, x);
  assert(x0 == 3.0);

  mcmclib_mh_q_free(q);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_rng_free(r);

  return 0;
}
