/**Test multivariate normal mixture distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mvnorm.h>
#include <mixnorm.h>

#define TOL 1e-6

#define DIM 10
#define V0 2.0
#define MU0 1.5

static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

int main(int argc, char** argv) {
  static mcmclib_mvnorm_lpdf* pis[2];

  /*set a non-trivial covariance matrix*/
  gsl_matrix* Sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma);
  gsl_matrix_scale(Sigma, V0);
  for(int i=0; i<DIM; i++) for(int j=i+1; j < DIM; j++) {
      gsl_matrix_set(Sigma, i, j, 1.0);
      gsl_matrix_set(Sigma, j, i, 1.0);
  }

  gsl_vector* mu1 = gsl_vector_alloc(DIM);
  gsl_vector_set_all(mu1, MU0);
  pis[0] = mcmclib_mvnorm_lpdf_alloc(mu1, Sigma->data);
  gsl_vector* mu2 = gsl_vector_alloc(DIM);
  gsl_vector_set_all(mu2, MU0);
  pis[1] = mcmclib_mvnorm_lpdf_alloc(mu2, Sigma->data);

  gsl_vector* w = gsl_vector_alloc(2);
  gsl_vector_set_all(w, 0.5);
  mcmclib_mixnorm_lpdf* p = mcmclib_mixnorm_lpdf_alloc(w, pis);

  gsl_vector* x = gsl_vector_alloc(DIM);

  for(int i=0; i<20; i++) {
    gsl_vector_set_all(x, 6.0 * (double) i / 20.0 - 3.0);
    double lpdf = mcmclib_mixnorm_lpdf_compute(p, x);
    double lpdf_check = exp(mcmclib_mvnorm_lpdf_compute(pis[0], x)) * gsl_vector_get(w, 0);
    lpdf_check += exp(mcmclib_mvnorm_lpdf_compute(pis[1], x)) * gsl_vector_get(w, 1);
    lpdf_check = log(lpdf_check);
    assert(check_dequal(lpdf, lpdf_check));
  }

  gsl_vector_free(x);
  mcmclib_mixnorm_lpdf_free(p);
  mcmclib_mvnorm_lpdf_free(pis[1]);
  mcmclib_mvnorm_lpdf_free(pis[0]);
  gsl_vector_free(mu2);
  gsl_vector_free(mu1);
  gsl_matrix_free(Sigma);
}
