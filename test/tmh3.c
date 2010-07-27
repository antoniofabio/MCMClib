/**Test MH object on a discrete distribution, asymm. proposal*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mh.h>
#include "CuTest.h"

#define N 10000
#define PI1 0.8

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

#define TOL 1e-6

static double dtarget(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  if(x0 == 0.0)
    return log(1.0 - PI1);
  else if(x0 == 1.0)
    return log(PI1);
  else
    return log(0.0);
}

static double qd(mcmclib_mh_q* q, gsl_vector* x, gsl_vector* y) {
  q = NULL; x = NULL; /*keep compiler quiet*/
  if(gsl_vector_get(y,0) == 0.0)
    return log(1.0 / 3.0);
  return log(2.0 / 3.0);
}

static void sampler(mcmclib_mh_q* q, gsl_vector* x) {
  gsl_rng* r = q->r;
  double coin = gsl_rng_uniform(r);
  if(coin < 1.0/3.0)
    gsl_vector_set(x, 0, 0.0);
  else
    gsl_vector_set(x, 0, 1.0);
}

void Testmh3(CuTest* tc) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, 0.0);
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, sampler, qd, NULL, NULL);
  mcmclib_mh* s = mcmclib_mh_alloc(r, dtarget, NULL, q, x);

  double pi1 = 0.0;

  for(int n=0; n<N; n++) {
    mcmclib_mh_update(s);
    pi1 += x0;
  }
  pi1 /= (double) N;
  CuAssertDblEquals(tc, 0.795100, pi1, TOL);

  gsl_vector_free(x);
  gsl_rng_free(r);
  mcmclib_mh_free(s);
}
