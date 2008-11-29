/**Test OLEM-RAPT algorithm on a mixture target*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mixnorm.h>
#include <olem_rapt.h>

static const double beta = 0.5;
static const double V[] = {1.0, 1.0};
static const double MU[] = {-0.5, 0.5};
static const double rho[] = {-0.6, -0.6};
static const int DIM = 2;

#define N 10000
#define DIM 1
#define K 2
#define T0 100

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

int main(int argc, char** argv) {
  /*setup target distrib.*/
  double w[] = {beta, 1-beta};
  gsl_vector_view wv = gsl_vector_view_array(w, K);
  gsl_matrix* Sigma[K];
  gsl_vector* mu[K];
  mcmclib_mvnorm_lpdf* pik[K];
  for(int k=0; k<K; k++) {
    Sigma[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma[k]);
    for(int i=0; i<DIM; i++) for(int j=0; j<i; j++) {
	gsl_matrix_set(Sigma[k], i, j, rho[k]);
	gsl_matrix_set(Sigma[k], j, i, rho[k]);
      }
    gsl_matrix_scale(Sigma[k], V[k]);
    mu[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu[k], MU[k]);
    pik[k] = mcmclib_mvnorm_lpdf_alloc(mu[k], Sigma[k]->data);
  }
  mcmclib_mixnorm_lpdf* pi = mcmclib_mixnorm_lpdf_alloc(&(wv.vector), pik);

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

  gsl_vector* mu_hat[K];
  gsl_matrix* Sigma_hat[K];
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], MU[k] * 0.5);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
    gsl_matrix_scale(Sigma_hat[k], 0.5);
  }
  gsl_vector* w_hat = gsl_vector_alloc(K);
  gsl_vector_set_all(w_hat, 1.0 / (double) K);

  mcmclib_olemrapt* sampler = mcmclib_olemrapt_alloc(rng,
						 mcmclib_mixnorm_lpdf_compute, pi,
						 x, T0, sigma_whole,
						 w_hat, mu_hat, Sigma_hat);

  /*Main MCMC loop*/
  double sum_x = 0.0;
  double sum_x2 = 0.0;
  for(int n=0; n<N; n++) {
    mcmclib_olemrapt_update(sampler);
    mcmclib_olemrapt_update_proposals(sampler);
    sum_x += gsl_vector_get(x, 0);
    sum_x2 += pow(gsl_vector_get(x, 0), 2);
  }
  /*check results*/
  printf("%f\t%f\n", sum_x, sum_x2);
  assert(check_dequal(sum_x, -16.94049));
  assert(check_dequal(sum_x2, 13348.093474));
  mcmclib_mixem_online* em = sampler->em;
  printf("mu1 = %f\t, mu2=%f\n", em->mu[0]->data[0], em->mu[1]->data[0]);
  assert(check_dequal(em->mu[0]->data[0], -0.551998));
  assert(check_dequal(em->mu[1]->data[0], 0.536395));

  /*free memory*/
  for(int k=0; k<K; k++)
    gsl_matrix_free(sigma_local[k]);
  gsl_matrix_free(sigma_whole);
  gsl_vector_free(x);
  mcmclib_olemrapt_free(sampler);
  for(int k=0; k<K; k++) {
    gsl_vector_free(mu[k]);
    gsl_matrix_free(Sigma[k]);
    mcmclib_mvnorm_lpdf_free(pik[k]);
  }
  mcmclib_mixnorm_lpdf_free(pi);

  return 0;
}
