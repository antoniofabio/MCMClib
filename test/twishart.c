/**Test Wishart distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lpdf_wishart.h>

#define TOL 1e-6

#define DIM 10
#define V0 2.0
#define P 2

static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

mcmclib_wishart_lpdf* p;
gsl_vector* x;

double lpdf(double s) {
  gsl_vector_scale(x, s);
  double ans = mcmclib_wishart_lpdf_compute(p, x);
  gsl_vector_scale(x, 1.0 / s);
  return ans;
}

int main(int argc, char** argv) {
  /*set a non-trivial location matrix*/
  gsl_matrix* V = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(V);
  gsl_matrix_scale(V, V0);
  for(int i=0; i<DIM; i++) for(int j=i+1; j < DIM; j++) {
      gsl_matrix_set(V, i, j, 1.0);
      gsl_matrix_set(V, j, i, 1.0);
  }
  
  p = mcmclib_wishart_lpdf_alloc(V, P);
  x = gsl_vector_alloc(DIM * DIM);
  gsl_matrix_view X_v = gsl_matrix_view_vector(x, DIM, DIM);
  gsl_matrix* X = &(X_v.matrix);
  gsl_matrix_set_all(X, 0.2);
  for(int i=0; i<DIM; i++)
    gsl_matrix_set(X, i, i, 1.0);

  assert(check_dequal(lpdf(1.0), -7.152627));
  assert(check_dequal(lpdf(0.5), -29.549142));
  assert(check_dequal(lpdf(2.0), 13.380252));

  mcmclib_wishart_lpdf_free(p);
  gsl_vector_free(x);
  gsl_matrix_free(V);
}
