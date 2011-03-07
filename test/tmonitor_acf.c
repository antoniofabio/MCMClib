#include <gsl/gsl_matrix.h>
#include <monitor.h>
#include "CuTest.h"

#define TOL 1e-6
#define MAX_LAG 1

#define x0 (x->data[0])

static gsl_vector* x;
static gsl_vector* y;
static mcmclib_monitor_acf_h m;
static gsl_matrix* ACF;

static void append(double xi) {
  gsl_vector_set(x, 0, xi);
  mcmclib_monitor_acf_update(m, x);
}

static double acf(size_t l) {
  mcmclib_monitor_acf_get(m, ACF);
  return gsl_matrix_get(ACF, l, 0);
}

void Testmonitor_acf(CuTest* tc) {
  x = gsl_vector_alloc(1);
  y = gsl_vector_alloc(1);
  m = mcmclib_monitor_acf_alloc(1, MAX_LAG);
  ACF = gsl_matrix_alloc(MAX_LAG+1, 1);

  append(1.0);
  append(-1.0);
  CuAssertDblEquals(tc, 1.0, acf(0), TOL);

  append(1.0);
  double mean = 1.0/3.0;
  double var = 1.0 - mean*mean;
  CuAssertDblEquals(tc, var, acf(0), TOL);

  double acf_check = -2.0/3.0;
  CuAssertDblEquals(tc, acf_check, acf(1), TOL);

  append(-1.0);
  mean = 0.0;
  var = 1.0;
  CuAssertDblEquals(tc, 1.0, acf(0), TOL);
  CuAssertDblEquals(tc, -3.0/4.0, acf(1), TOL);

  gsl_vector* iact = gsl_vector_alloc(1);
  mcmclib_monitor_acf_get(m, ACF);
  mcmclib_iact_from_acf(ACF, iact);
  CuAssertDblEquals(tc, -0.5, gsl_vector_get(iact, 0), TOL);

  gsl_vector_free(iact);
  gsl_matrix_free(ACF);
  mcmclib_monitor_acf_free(m);
  gsl_vector_free(x);
  gsl_vector_free(y);
}
