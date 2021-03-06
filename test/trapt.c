/**Test RAPT algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "rapt.c"
#include "CuTest.h"

#define N 1000
#define DIM 1
#define K 2
/*burn-in*/
#define T0 100

#define TOL 1e-6

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static size_t which_region(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  return x0 < 0.5 ? 0 : 1;
}

static double dunif(void* ignore, const gsl_vector* x) {
  ignore = NULL; /*keep compiler quiet*/
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

void Testrapt(CuTest* tc) {
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_matrix* sigma_whole = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma_whole);
  gsl_matrix* sigma_local[K];
  for(int k=0; k<K; k++) {
    sigma_local[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(sigma_local[k]);
  }

  double means[K];
  double variances[K];
  double nk[K];
  for(int k=0; k<K; k++) {
    means[k] = 0.0;
    variances[k] = 0.0;
    nk[k] = 0.0;
  }
  double mean = 0.0;
  double variance = 0.0;

  mcmclib_amh* s = mcmclib_rapt_alloc(rng,
				      dunif, NULL, /*target distrib.*/
				      x, T0,
				      sigma_whole, K, sigma_local,
				      which_region, NULL, NULL);
  rapt_suff* suff = (rapt_suff*) s->suff;

  /*Main MCMC loop*/
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  gsl_vector* which_region_n = gsl_vector_alloc(N);
  for(size_t n=0; n<N; n++) {
    mcmclib_amh_update(s);

    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x);
    gsl_vector_set(which_region_n, n, (double) which_region(NULL, x));
    means[which_region(NULL, x)] += x0;
    variances[which_region(NULL, x)] += x0 * x0;
    nk[which_region(NULL, x)] += 1.0;
    mean += x0;
    variance += x0 * x0;
  }

  /*compute means and variances*/
  mean /= (double) N;
  variance = variance / ((double) N) - (mean * mean);
  for(size_t k=0; k<K; k++) {
    means[k] /= nk[k];
    variances[k] = (variances[k] / nk[k]) - (means[k] * means[k]);
  }

  /*check results*/
  CuAssertDblEquals(tc, mean, v0(suff->global_mean), TOL);
  CuAssertDblEquals(tc, variance, m00(suff->global_variance), TOL);
  static char kmsg[3];
  for(size_t k=0; k<K; k++) {
    sprintf(kmsg, "%zd", k);
    CuAssertDblEquals_Msg(tc, kmsg, nk[k], gsl_vector_get(suff->n, k), TOL);
    CuAssertDblEquals_Msg(tc, kmsg, means[k], v0(suff->means[k]), TOL);
    CuAssertDblEquals_Msg(tc, kmsg, variances[k], m00(suff->variances[k]), TOL);
  }

  /*free memory*/
  gsl_matrix_free(X);
  for(int k=0; k<K; k++)
    gsl_matrix_free(sigma_local[k]);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_amh_free(s);
  gsl_rng_free(rng);
  gsl_vector_free(which_region_n);
}
