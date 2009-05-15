/**Test linear regression model */
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <lm.h>

#define K 3
#define N 10

int main(int argc, char** argv) {
  gsl_vector *y, *beta, *sigmasq;
  gsl_matrix *X;
  X = gsl_matrix_alloc(N, K);
  gsl_matrix_set_zero(X);
  for(int i=0; i<K; i++)
    gsl_matrix_set(X, i, i, 1.0);
  y = gsl_vector_alloc(N);
  gsl_vector_set_all(y, 1.0);
  beta = gsl_vector_alloc(K);
  sigmasq = gsl_vector_alloc(1);
  mcmclib_lm* lm = mcmclib_lm_alloc(X, y, beta, sigmasq);
  gsl_matrix_free(X);

  mcmclib_lm_free(lm);
  gsl_vector_free(y);
  gsl_vector_free(beta);
  gsl_vector_free(sigmasq);
}
