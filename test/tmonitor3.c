/**Test monitor object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <monitor.h>
#include "CuTest.h"

#define TOL 1e-6

void Testmonitor3(CuTest* tc) {
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, -1.0);
  mcmclib_monitor_h p = mcmclib_monitor_alloc(x);
  mcmclib_monitor_update(p);
  gsl_vector* tmp = gsl_vector_alloc(x->size);

  /*check that acceptance rate never goes above 1...*/
  for(int i = 0; i < 1000; i++) {
    gsl_vector_set(x, 0, (double) i);
    mcmclib_monitor_update(p);
    mcmclib_monitor_get_ar(p, tmp);
    double ar = gsl_vector_get(tmp, 0);
    CuAssertTrue(tc, gsl_finite(ar));
    CuAssertTrue(tc, ar <= 1.0);
    CuAssertDblEquals(tc, 1.0, ar, TOL);
  }

  mcmclib_monitor_free(p);
  gsl_vector_free(tmp);
  gsl_vector_free(x);
}
