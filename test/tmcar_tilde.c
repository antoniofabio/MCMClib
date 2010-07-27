/**Test spatial normal distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mcar_tilde.h>
#include <matrix.h>
#include "CuTest.h"

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static int is_pos_def(gsl_matrix* A) {
  if(A->size1 != A->size2)
    return 0;
  size_t n = A->size1;
  for(size_t i=0; i<(n-1); i++) {
    if(gsl_matrix_get(A, i, i) <= 0.0)
      return 0;
    for(size_t j=(i+1); j<n; j++)
      if(fabs(gsl_matrix_get(A, i, j) - gsl_matrix_get(A, j, i)) >= TOL)
	return 0;
  }
  return 1;
}

#define P 2
#define N 3

static mcmclib_mcar_tilde_lpdf* p;

/* declare as near regions 'i' and 'j' */
#define DECL_AD(i, j) if(1) {			\
    gsl_matrix_set(W, i, j, 1.0);		\
    gsl_matrix_set(W, j, i, 1.0);		\
  }

static gsl_vector* x;
static double lpdf(double s) {
  gsl_vector_set_all(x, s);
  return mcmclib_mcar_tilde_lpdf_compute(p, x);
}

void Testmcar_tilde(CuTest* tc) {
  gsl_vector* mu = gsl_vector_alloc(N * P);
  gsl_vector_set_zero(mu);
  gsl_matrix* W = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(W);
  for(size_t i=0; i<(N-1); i++)
    DECL_AD(i, i+1);
  
  p = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  for(size_t i=0; i<N; i++) {
    int count = 0;
    for(size_t j=0; j<N; j++)
      count += (gsl_matrix_get(W, i, j) == 1.0) ? 1 : 0;
    CuAssertTrue(tc, gsl_vector_get(p->m, i) == (double) count);
  }

  gsl_vector_set_all(p->alphasigmag, 0.0);
  gsl_vector_set_all(p->alpha12sigma, -1.0);

  mcmclib_mcar_tilde_lpdf_update_vcov(p);

  gsl_vector_set_all(p->alpha12sigma, -0.7);
  mcmclib_mcar_tilde_lpdf_update_blocks(p);

  x = gsl_vector_alloc(N*P);
  gsl_vector_set_all(p->alpha12sigma, -2.0);
  double l1 = lpdf(0.0);
  CuAssertTrue(tc, gsl_finite(l1));
  CuAssertTrue(tc, gsl_finite(lpdf(1.0)));
  CuAssertTrue(tc, l1 == lpdf(0.0)); /*check for side-effects*/
  gsl_vector_set_all(p->alpha12sigma, 0.1);
  CuAssertTrue(tc, !gsl_finite(lpdf(1.0)));

  gsl_vector_free(x);

  mcmclib_mcar_tilde_lpdf_free(p);
  gsl_matrix_free(W);
  gsl_vector_free(mu);
}
