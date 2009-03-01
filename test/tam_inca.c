/**Test AM-INCA algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>
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
  for(int n=0; n<N; n++) {
    mcmclib_am_inca_update(s);

    for(int m=0; m<M; m++) {
      mean += v0(x[m]);
      variance += v0(x[m]) * v0(x[m]);
    }
  }


  /*check results*/
  assert(s->inca->amh->n == (N * M));
  mcmclib_gauss_am_suff* suff = (mcmclib_gauss_am_suff*) s->inca->amh->suff;

  assert(check_dequal(mean, v0(suff->sum_x)));
  assert(check_dequal(variance, m00(suff->sum_xx)));
  double eps = m00(suff->Sigma_eps);
  mcmclib_gauss_mrw_gamma* gamma = (mcmclib_gauss_mrw_gamma*) s->inca->amh->mh->q->gamma;
  mean /= (double) (N * M);
  variance = variance / ((double) (N * M)) - (mean * mean);
  assert(check_dequal(fix(variance, eps), m00(gamma->Sigma)));

  /*free memory*/
  gsl_matrix_free(sigma);
  for(int m=0; m<M; m++)
    gsl_vector_free(x[m]);
  mcmclib_am_inca_free(s);

  return 0;
}
