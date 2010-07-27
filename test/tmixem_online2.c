/**Test online EM mixture fitting algorithm*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mvnorm.h>
#include <mixem_online.h>
#include "CuTest.h"

static const size_t DIM = 2;
static const size_t N = 5;
static const size_t K = 2;

static size_t count_errors = 0;
static void test_handler (const char * reason,
			  const char * file,
			  int line,
			  int gsl_errno) {
  count_errors++;
}

void Testmixem_online2(CuTest* tc) {
  /*init data*/
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  gsl_matrix_set_all(X, 1.0);

  /*fit mixture by online EM*/
  gsl_vector* mu_hat[K];
  gsl_matrix* Sigma_hat[K];
  for(size_t k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], 1.0);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_zero(Sigma_hat[k]);
  }
  gsl_vector* w_hat = gsl_vector_alloc(K);
  gsl_vector_set_all(w_hat, 1.0 / (double) K);

  mcmclib_mixem_online* m = mcmclib_mixem_online_alloc(mu_hat,
						       Sigma_hat,
						       w_hat,
						       0.5, 2);
  gsl_set_error_handler(test_handler);
  for(size_t n=0; n<N; n++) {
    gsl_vector_view rv = gsl_matrix_row(X, n);
    mcmclib_mixem_online_update(m, &(rv.vector));
  }
  CuAssertTrue(tc, count_errors == N);
  mcmclib_mixem_online_free(m);
  gsl_matrix_free(X);

  /*free memory*/
  for(size_t k=0; k<K; k++) {
    gsl_matrix_free(Sigma_hat[k]);
    gsl_vector_free(mu_hat[k]);
  }
  gsl_vector_free(w_hat);
}
