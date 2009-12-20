/**Test AMH object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <amh.h>

#define N 10000
#define PI1 0.8

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

static double dtarget(void* ignore, gsl_vector* x) {
  if(x0 == 0.0)
    return log(1.0 - PI1);
  else if(x0 == 1.0)
    return log(PI1);
  else
    return log(0.0);
}

static double qd(mcmclib_mh_q* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

static void sampler(mcmclib_mh_q* q, gsl_vector* x) {
  gsl_rng* r = q->r;
  double coin = gsl_rng_uniform(r);
  if(coin < 0.5)
    gsl_vector_set(x, 0, 0.0);
  else
    gsl_vector_set(x, 0, 1.0);
}

static void update_gamma(mcmclib_amh* p) {
  //DOES NOTHING
}

int main(int argc, char** argv) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, -0.5);
  mcmclib_amh* s = mcmclib_amh_alloc(mcmclib_mh_alloc(r, dtarget, NULL,
						      mcmclib_mh_q_alloc(r, sampler, qd, NULL, NULL), x),
				     &gamma, update_gamma, NULL);

  double pi1 = 0.0;
  for(int n=0; n<N; n++) {
    mcmclib_amh_update(s);
    pi1 += x0;
  }
  pi1 /= (double) N;
  assert(0.795900 == pi1);

  gsl_vector_free(x);
  gsl_rng_free(r);
  mcmclib_amh_free(s);

  return 0;
}
