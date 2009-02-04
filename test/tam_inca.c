/**Test base AM-INCA algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <am_inca.h>

#define N 1000
#define DIM 1
#define M 3 /*number of parallel chains*/
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

  gsl_vector* x[M];
  for(int m=0; m<M; m++) {
    x[m] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(x[m], gsl_rng_uniform(rng));
  }

  gsl_matrix* sigma = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(sigma);

  double mean = 0.0;
  double variance = 0.0;

  mcmclib_am_inca* s = mcmclib_am_inca_alloc(rng,
					     dunif, NULL, /*target distrib.*/
					     x, M, T0, sigma);

  /*Main MCMC loop*/
  gsl_matrix* X = gsl_matrix_alloc(N * M, DIM);
  for(int n=0; n<N; n++) {
    mcmclib_am_inca_update(s);

    for(int m=0; m<M; m++) {
      int n1 = n * M + m;
      gsl_vector_view Xn = gsl_matrix_row(X, n1);
      gsl_vector_memcpy(&(Xn.vector), x[m]);
      mean += v0(x[m]);
      variance += v0(x[m]) * v0(x[m]);
    }
  }

  /*compute mean and variance*/
  mean /= (double) (N * M);
  variance = variance / ((double) (N * M)) - (mean * mean);

  /*check results*/
  assert(s->t == (N * M));
  assert(check_dequal(mean, v0(s->global_mean)));
  assert(check_dequal(variance, m00(s->global_variance)));
  double eps = m00(s->Sigma_eps);
  assert(check_dequal(fix(variance, eps), m00(s->sm[M-1]->sigma_prop)));

  /*free memory*/
  gsl_matrix_free(sigma);
  for(int m=0; m<M; m++)
    gsl_vector_free(x[m]);
  mcmclib_am_inca_free(s);

  return 0;
}
