/**Test RAPT proposal kernel*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt_q.h>
#include "CuTest.h"

#define DIM 1
#define K 2

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

static size_t which_region(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  return x0 < 0.5 ? 0 : 1;
}

#define TOL 1e-6

void Testrapt_q(CuTest* tc) {
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector* y = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);
  gsl_matrix* sigma_local[K];
  for(size_t k=0; k<K; k++) {
    sigma_local[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(sigma_local[k]);
  }

  mcmclib_mh_q* q = mcmclib_rapt_q_alloc(rng,
					 sigma_whole, K, sigma_local,
					 which_region, NULL, NULL);

  gsl_vector_memcpy(y, x);
  mcmclib_mh_q_sample(q, x);
  CuAssertDblEquals(tc, -1.048402, mcmclib_mh_q_logd(q, x, y), TOL);

  /*free memory*/
  for(size_t k=0; k<K; k++)
    gsl_matrix_free(sigma_local[k]);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  gsl_vector_free(y);
  mcmclib_mh_q_free(q);
  gsl_rng_free(rng);
}
