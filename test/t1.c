/**Test multivariate normal distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mvnorm.h>

#define TOL 1e-6

#define DIM 10
#define V0 2.0
#define MU0 1.5

static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

int main(int argc, char** argv) {
  /*set a non-trivial covariance matrix*/
  gsl_matrix* Sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma);
  gsl_matrix_scale(Sigma, V0);
  for(int i=0; i<DIM; i++) for(int j=i+1; j < DIM; j++) {
      gsl_matrix_set(Sigma, i, j, 1.0);
      gsl_matrix_set(Sigma, j, i, 1.0);
  }

  gsl_vector* mu = gsl_vector_alloc(DIM);
  gsl_vector_set_all(mu, MU0);

  mcmclib_mvnorm_lpdf* p = mcmclib_mvnorm_lpdf_alloc(mu, Sigma->data);
  mcmclib_mvnorm_lpdf* p_nochol = mcmclib_mvnorm_lpdf_alloc(mu, Sigma->data);
  mcmclib_mvnorm_lpdf_chol(p_nochol);
  mcmclib_mvnorm_lpdf* p_noinv = mcmclib_mvnorm_lpdf_alloc(mu, Sigma->data);
  mcmclib_mvnorm_lpdf_inverse(p_noinv);
  gsl_vector* x = gsl_vector_alloc(DIM);

  gsl_matrix* check_m = gsl_matrix_alloc(20, 2);
  FILE* check_f = fopen("t1.dat.check", "r");
  gsl_matrix_fscanf(check_f, check_m);
  fclose(check_f);

  for(int i=0; i<20; i++) {
    gsl_vector_set_all(x, gsl_matrix_get(check_m, i, 0));
    double lpdf = mcmclib_mvnorm_lpdf_compute(p, x);
    double lpdf_check = gsl_matrix_get(check_m, i, 1);
    assert(check_dequal(lpdf, lpdf_check));
    assert(check_dequal(mcmclib_mvnorm_lpdf_compute_nochol(p_nochol, x), lpdf_check));
    assert(check_dequal(mcmclib_mvnorm_lpdf_compute_noinv(p_noinv, x), lpdf_check));
  }

  gsl_matrix_free(check_m);
  mcmclib_mvnorm_lpdf_free(p_noinv);
  mcmclib_mvnorm_lpdf_free(p_nochol);
  mcmclib_mvnorm_lpdf_free(p);

  /*negative correlations*/
  for(int i=0; i<DIM; i++) {
    gsl_matrix_set(Sigma, i, i, 1.0);
    for(int j=0; j<i; j++){
      gsl_matrix_set(Sigma, i, j, -1.0 / (double)DIM);
      gsl_matrix_set(Sigma, j, i, -1.0 / (double)DIM);
    }
  }
  p = mcmclib_mvnorm_lpdf_alloc(mu, Sigma->data);
  gsl_vector_set_all(x, -1.0);
  assert(check_dequal(mcmclib_mvnorm_lpdf_compute(p, x), -320.966989));
  for(int i=0; i<DIM; i++)
    gsl_vector_set(x, i, pow(-1, i+1));
  assert(check_dequal(mcmclib_mvnorm_lpdf_compute(p, x), -125.512443));
  gsl_vector_free(x);
  gsl_vector_free(mu);
  gsl_matrix_free(Sigma);
}
