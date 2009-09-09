/**Test GAUSS-AM algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>

#define N 1000
#define DIM 1
/*burn-in*/
#define T0 100
#define SF (2.38*2.38/(double) DIM)

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

static double dunif(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  return log(0.0);
}

static double fix(double in, double correction) {
  return (in + correction) * SF;
}

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, gsl_rng_uniform(rng));

  gsl_matrix* sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma);

  double sum_x = 0.0;
  double sum_xx = 0.0;

  mcmclib_amh* s = mcmclib_gauss_am_alloc(rng,
					  dunif, NULL, /*target distrib.*/
					  x, sigma, T0);

  /*Main MCMC loop*/
  for(int n=0; n<N; n++) {
    mcmclib_amh_update(s);

    sum_x += x0;
    sum_xx += x0 * x0;
  }

  /*check results*/
  assert(s->n == N);
  mcmclib_gauss_am_suff* suff = (mcmclib_gauss_am_suff*) s->suff;

  assert(check_dequal(sum_x, v0(suff->sum_x)));
  assert(check_dequal(sum_xx, m00(suff->sum_xx)));
  double eps = m00(suff->Sigma_eps);
  mcmclib_gauss_mrw_gamma* gamma = (mcmclib_gauss_mrw_gamma*) s->mh->q->gamma;
  double mean = sum_x / (double) N;
  double variance = sum_xx / (double) N - mean * mean;
  assert(check_dequal(fix(variance, eps), m00(gamma->Sigma)));

  /*check against uniform distribution*/
  assert(check_dequal(0.461706, mean));
  assert(check_dequal(0.076635, variance));

  /*free memory*/
  gsl_matrix_free(sigma);
  gsl_vector_free(x);
  mcmclib_gauss_am_free(s);
  gsl_rng_free(rng);

  return 0;
}
