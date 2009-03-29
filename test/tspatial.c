/**Test spatial normal distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <spatial.h>

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define N 9
#define MU 1.5
#define RHO 2.0
#define SIGMA2 2.0
#define TAUSQ 0.5

gsl_vector* x;
mcmclib_spatial_lpdf* p;

double lpdf(double s) {
  gsl_vector_set_all(x, s);
  return mcmclib_spatial_lpdf_compute(p, x);
}

int main(int argc, char** argv) {
  gsl_vector* mu = gsl_vector_alloc(N);
  gsl_vector_set_all(mu, MU);
  gsl_vector* rho = gsl_vector_alloc(1);
  gsl_vector_set(rho, 0, RHO);
  gsl_vector* sigma = gsl_vector_alloc(1);
  gsl_vector_set(sigma, 0, SIGMA2);
  gsl_vector* tausq = gsl_vector_alloc(1);
  gsl_vector_set(tausq, 0, TAUSQ);
  gsl_matrix* D = gsl_matrix_alloc(N, N);

  gsl_matrix* xy = gsl_matrix_alloc(N, 2);
  int k=0;
  for(int i=0; i<3; i++) for(int j=0; j<3; j++){
      gsl_matrix_set(xy, k, 0, (double) i);
      gsl_matrix_set(xy, k, 1, (double) j);
      k++;
    }
  mcmclib_spatial_distances(D, xy);
  gsl_matrix_free(xy);

  p = mcmclib_spatial_lpdf_alloc(mu, rho, sigma, tausq, D);
  x = gsl_vector_alloc(N);
  assert(check_dequal(lpdf(1.0), -10.253884));
  assert(check_dequal(lpdf(0.0), -11.491107));
  assert(check_dequal(lpdf(-1.0),-13.965553));
  assert(check_dequal(lpdf(1.5), -10.099231));
  gsl_vector_free(x);

  mcmclib_spatial_lpdf_free(p);
  gsl_matrix_free(D);
  gsl_vector_free(mu);
  gsl_vector_free(rho);
  gsl_vector_free(sigma);
  gsl_vector_free(tausq);
}
