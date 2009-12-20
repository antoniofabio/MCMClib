/**Test RAPT proposal kernel*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt_q.h>

#define DIM 1
#define K 2

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

static int which_region(void* ignore, gsl_vector* x) {
  return x0 < 0.5 ? 0 : 1;
}

static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

int main(int argc, char** argv) {
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector* y = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);
  gsl_matrix* sigma_local[K];
  for(int k=0; k<K; k++) {
    sigma_local[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(sigma_local[k]);
  }

  mcmclib_mh_q* q = mcmclib_rapt_q_alloc(rng,
					 dunif, NULL, /*target distrib.*/
					 sigma_whole, K, sigma_local,
					 which_region, NULL);

  gsl_vector_memcpy(y, x);
  mcmclib_mh_q_sample(q, x);
  assert(check_dequal(mcmclib_mh_q_logd(q, x, y), -1.048402));

  /*free memory*/
  for(int k=0; k<K; k++)
    gsl_matrix_free(sigma_local[k]);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  gsl_vector_free(y);
  mcmclib_mh_q_free(q);
  gsl_rng_free(rng);

  return 0;
}
