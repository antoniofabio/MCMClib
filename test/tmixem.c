/**Test mixture fitting em algorithm*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mvnorm.h>
#include <mixem.h>

static const double beta = 0.8;
static const double V[] = {1.0, 4.0};
static const double MU[] = {-3.0, 3.0};
static const int DIM = 3;

#define N 1000

#define K 2

static int sample(gsl_rng* r, gsl_vector* probs);

int main(int argc, char** argv) {
  /*setup target distrib. parameters*/
  double w[] = {beta, 1-beta};
  gsl_vector_view wv = gsl_vector_view_array(w, K);
  gsl_matrix* Sigma[K];
  gsl_vector* mu[K];
  for(int k=0; k<K; k++) {
    Sigma[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma[k]);
    gsl_matrix_scale(Sigma[k], V[k]);
    mu[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu[k], MU[k]);
  }

  /*generate random data*/
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  for(int n=0; n<N; n++) {
    gsl_vector_view rv = gsl_matrix_row(X, n);
    gsl_vector* r = &(rv.vector);
    int k = sample(rng, &(wv.vector));
    mcmclib_mvnorm(rng, Sigma[k], r);
    gsl_vector_add(r, mu[k]);
  }
  FILE* out_X = fopen("tmixem_X.dat", "w");
  gsl_matrix_fprintf(out_X, X, "%f");
  fclose(out_X);

  /*fit mixture*/
  gsl_vector* mu_hat[K];
  gsl_matrix* Sigma_hat[K];
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], MU[k]);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
  }
  gsl_matrix* P = gsl_matrix_alloc(N, K);
  gsl_vector* w_hat = gsl_vector_alloc(K);
  mcmclib_mixem_fit(X, K, mu_hat, Sigma_hat, P, w_hat, 0.001);

  /*print out estimation results*/
  gsl_vector_fprintf(stdout, w_hat, "%f");
  FILE* out_mu = fopen("tmixem_mu.dat", "w");
  FILE* out_Sigma = fopen("tmixem_Sigma.dat", "w");
  for(int k=0; k<K; k++) {
    gsl_vector_fprintf(out_mu, mu_hat[k], "%f");
    gsl_matrix_fprintf(out_Sigma, Sigma_hat[k], "%f");
  }
  fclose(out_Sigma);
  fclose(out_mu);

  return 0;
}

static int sample(gsl_rng* r, gsl_vector* probs) {
  int size = probs->size;
  double cum_sum = 0.0;
  double who = gsl_rng_uniform(r);
  for(int which=0; which<size; which++) {
    if(who < (cum_sum += gsl_vector_get(probs, which)))
      return(which);
  }
  return(size-1);
}
