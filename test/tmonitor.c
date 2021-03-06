/**Test monitor object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <monitor.h>
#include "CuTest.h"

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

void Testmonitor(CuTest* tc) {
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, 2.0);
  mcmclib_monitor_h p = mcmclib_monitor_alloc(x);
  gsl_vector* tmp = gsl_vector_alloc(1);

  /*check state after 1st update*/
  mcmclib_monitor_update(p);
  mcmclib_monitor_get_means(p, tmp);
  CuAssertDblEquals(tc, x0, v0(tmp), TOL);
  mcmclib_monitor_get_vars(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_ar(p, tmp);
  CuAssertTrue(tc, !gsl_finite(v0(tmp)));
  mcmclib_monitor_get_msjd(p, tmp);
  CuAssertTrue(tc, !gsl_finite(v0(tmp)));

  /*check state after 1 rejection*/
  mcmclib_monitor_update(p);

  mcmclib_monitor_get_means(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), x0));
  mcmclib_monitor_get_vars(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_ar(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_msjd(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 0.0));

  /*check state after 1 rejection and 1 acceptance*/
  gsl_vector_set(x, 0, -4.0);
  mcmclib_monitor_update(p);

  mcmclib_monitor_get_means(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_vars(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 8.0));
  mcmclib_monitor_get_ar(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 0.5));
  mcmclib_monitor_get_msjd(p, tmp);
  CuAssertTrue(tc, check_dequal(v0(tmp), 18.0)); /*(2+4)^2 / 2*/

  mcmclib_monitor_free(p);
  gsl_vector_free(tmp);
  gsl_vector_free(x);
}
