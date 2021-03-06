/**Test AMH object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <amh.h>
#include "CuTest.h"

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

static double dtarget(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  if((x0 < -10.0))
    return log(0.5);
  return log(0.0);
}

static void sampler(mcmclib_mh_q* p, gsl_vector* x) {
  double* o = (double*) p->gamma;
  gsl_vector_set(x, 0, x0 + (*o));
}

static void update_suff(void* p, const gsl_vector* x) {
  int *suff = (int*) p;
  suff[0]++;
}

void Testamh(CuTest *tc) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, 0.5);
  double inc = 0.2;
  int suff = 0;
  mcmclib_amh* s = mcmclib_amh_alloc(mcmclib_mh_alloc(r, dtarget, NULL,
						      mcmclib_mh_q_alloc(r, sampler, NULL, &inc, NULL), x), 0,
				     &suff, NULL, update_suff, NULL);
  for(int n=0; n<10; n++) {
    mcmclib_amh_update(s);
  }
  CuAssertIntEquals(tc, 10, suff);

  gsl_vector_free(x);
  gsl_rng_free(r);
  mcmclib_amh_free(s);
}
