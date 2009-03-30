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

static double qd(void* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

static void sampler(void* state, gsl_vector* x) {
  mcmclib_mh_q* q = (mcmclib_mh_q*) state;
  gsl_rng* r = (gsl_rng*) q->gamma;
  double coin = gsl_rng_uniform(r);
  if(coin < 0.5)
    gsl_vector_set(x, 0, 0.0);
  else
    gsl_vector_set(x, 0, 1.0);
}

static void update_gamma(void* p) {
  //DOES NOTHING
}

int main(int argc, char** argv) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, -0.5);
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, sampler, r, qd, NULL, r);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, dtarget, NULL, q, x);
  mcmclib_amh* s = mcmclib_amh_alloc(mh, NULL, update_gamma);

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
  mcmclib_mh_free(mh);
  mcmclib_mh_q_free(q);

  return 0;
}
