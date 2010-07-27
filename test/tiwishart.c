/**Test Inverse-Wishart distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lpdf_iwishart.h>
#include "CuTest.h"

#define TOL 1e-6

#define DIM 10
#define P DIM

static mcmclib_iwishart_lpdf* p;
static gsl_vector* x;

static double lpdf(double s) {
  gsl_vector* y = gsl_vector_alloc(DIM * DIM);
  gsl_vector_memcpy(y, x);
  gsl_vector_scale(y, s);
  double ans = mcmclib_iwishart_lpdf_compute(p, y);
  gsl_vector_free(y);
  return ans;
}

void Testiwishart(CuTest* tc) {
  /*set a non-trivial location matrix*/
  gsl_matrix* V = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(V);
  gsl_matrix_add_constant(V, 1.0);
  
  p = mcmclib_iwishart_lpdf_alloc(V, P);
  x = gsl_vector_alloc(DIM * DIM);
  gsl_matrix_view X_v = gsl_matrix_view_vector(x, DIM, DIM);
  gsl_matrix* X = &(X_v.matrix);
  gsl_matrix_set_identity(X);
  gsl_matrix_add_constant(X, 1.0);

  /*check for side-effects*/
  double tmp = lpdf(1.0);
  CuAssertTrue(tc, tmp == lpdf(1.0));

  CuAssertDblEquals(tc, -18.188424, lpdf(1.0), TOL);
  CuAssertDblEquals(tc, 49.59203, lpdf(0.5), TOL);
  CuAssertDblEquals(tc, -88.468878, lpdf(2.0), TOL);

  mcmclib_iwishart_lpdf_free(p);
  gsl_vector_free(x);
  gsl_matrix_free(V);
}
