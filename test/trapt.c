/**Test base RAPT algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>

#define N 1000
#define DIM 1
#define K 2
/*burn-in*/
#define T0 100

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static int which_region(gsl_vector* x, void* ignore) {
  return x0 < 0.5 ? 0 : 1;
}

static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

int main(int argc, char** argv) {

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

  mcmclib_rapt* s = mcmclib_rapt_alloc(rng,
				       dunif, NULL, /*target distrib.*/
				       x, T0,
				       sigma_whole, K, sigma_local,
				       which_region, NULL);

  /*Main MCMC loop*/
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  gsl_vector* which_region_n = gsl_vector_alloc(N);
  for(int n=0; n<N; n++) {
    mcmclib_rapt_update(s);
    mcmclib_rapt_update_proposals(s);

    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x);
    gsl_vector_set(which_region_n, n, s->which_region_x);
    means[s->which_region_x] += x0;
    variances[s->which_region_x] += x0 * x0;
    nk[s->which_region_x] += 1.0;
    mean += x0;
    variance += x0 * x0;
  }

  /*compute means and variances*/
  mean /= (double) N;
  variance = variance / ((double) N) - (mean * mean);
  for(int k=0; k<K; k++) {
    means[k] /= nk[k];
    variances[k] = (variances[k] / nk[k]) - (means[k] * means[k]);
  }

  /*check results*/
  assert(s->t == N);
  assert(s->t0 == T0);
  assert(check_dequal(mean, v0(s->global_mean)));
  assert(check_dequal(variance, m00(s->global_variance)));
  for(int k=0; k<K; k++) {
    assert(check_dequal(nk[k], gsl_vector_get(s->n, k)));
    assert(check_dequal(means[k], v0(s->means[k])));
    assert(check_dequal(variances[k], m00(s->variances[k])));
  }

  /*free memory*/
  for(int k=0; k<K; k++)
    gsl_matrix_free(sigma_local[k]);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_rapt_free(s);

  return 0;
}
