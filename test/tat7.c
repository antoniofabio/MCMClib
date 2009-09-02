/**Test AT7 algorithm on a dumb target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <at7.h>
#include <monitor.h>

static const double beta = 0.5;
static const double V[] = {1.0, 1.0};
static const double MU[] = {0.2, 0.8};

#define N 500
#define DIM 1
#define K 2
#define T0 250
#define SF (2.38*2.38/(double) DIM)

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)
#define m00(m) gsl_matrix_get(m, 0, 0)

#define TOL 1e-6
int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static double dunif(void* ignore, gsl_vector* x) {
  for(int i=0; i < x->size; i++) {
    double xi = gsl_vector_get(x, i);
    if((xi < 0.0) || (xi > 1.0))
      return log(0.0);
  }
  return log(1.0);
}

#define INIT_CHAIN							\
  if(1) {								\
  gsl_vector_set_all(x, 0.5);						\
  m = mcmclib_monitor_alloc(x);						\
  sampler = mcmclib_at7_alloc(rng,					\
			      dunif, NULL,				\
			      x, T0,					\
			      w_hat, mu_hat, Sigma_hat);		\
  }									\

#define RUN_CHAIN							\
  if(1) {								\
    for(int n=0; n<N; n++) {						\
      mcmclib_amh_update(sampler);					\
      mcmclib_monitor_update(m);					\
    }									\
    mcmclib_at7_free(sampler);						\
    mcmclib_monitor_free(m);						\
  }									\

#define TRY_DIM(dim)							\
  if(1) {								\
    x = gsl_vector_alloc(dim);						\
    for(int k=0; k<K; k++) {						\
      mu_hat[k] = gsl_vector_alloc(dim);				\
      gsl_vector_set_all(mu_hat[k], MU[k] * 0.5);			\
      Sigma_hat[k] = gsl_matrix_alloc(dim, dim);			\
      gsl_matrix_set_identity(Sigma_hat[k]);				\
    }									\
    gsl_vector_set_all(w_hat, 1.0 / (double) K);			\
    rng = gsl_rng_alloc(gsl_rng_default);				\
    INIT_CHAIN;								\
    RUN_CHAIN;								\
									\
    INIT_CHAIN;								\
    mcmclib_at7_set_sf_all(sampler, 0.2);				\
    RUN_CHAIN;								\
									\
    INIT_CHAIN;								\
    sf = gsl_vector_alloc(K);						\
    gsl_vector_set_all(sf, 2.0);					\
    mcmclib_at7_set_sf(sampler, sf);					\
    gsl_vector_free(sf);						\
    RUN_CHAIN;								\
  /*free memory*/							\
  gsl_vector_free(x);							\
  gsl_rng_free(rng);							\
  for(int k=0; k<K; k++) {						\
    gsl_vector_free(mu_hat[k]);						\
    gsl_matrix_free(Sigma_hat[k]);					\
  }									\
  }									\

int main(int argc, char** argv) {
  gsl_rng* rng;
  gsl_vector* x;
  mcmclib_amh* sampler;
  mcmclib_monitor* m;
  gsl_vector* sf;

  gsl_vector* mu_hat[K];
  gsl_matrix* Sigma_hat[K];
  gsl_vector* w_hat = gsl_vector_alloc(K);

  TRY_DIM(1);
  /*  TRY_DIM(2);
      TRY_DIM(3);*/

  gsl_vector_free(w_hat);
  return 0;
}
