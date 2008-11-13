/**RAPT example 4: mixture-based adaptive boundary*/
#include <stdio.h>
#include <memory.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>
#include <mixem.h>
#include <mvnorm.h>

/** Following arguments can be (almost) freely customized */
/*chain length*/
#define N 100000
/*state space dimension*/
#define DIM 2
/*absolute local mean value*/
#define MU0 1.5
/*local variances multipliers*/
static double V0[] = {1.0, 1.0};
/*burn in length*/
#define T0 ((DIM + DIM * (DIM-1) / 2) * 100)
/*update boundary every N0 iterations*/
#define N0 (N / 10)
/*scaling factor*/
#define SCALING_FACTOR (2.38 * 2.38 / (double) DIM)
/*starting local variance guess as scaling factor w.r.t. true value*/
#define SIGMA0_LOCAL 0.25
/*starting global variance guess*/
#define SIGMA0_GLOBAL pow(MU0 * 2, 2.0)
/***/

/*number of regions*/
#define K 2
/*boundary function and support data*/
static gsl_vector* mu_hat[K];
static gsl_matrix* Sigma_hat[K];
static mcmclib_mvnorm_lpdf* pi_hat[K];
int which_region(gsl_vector* x, void* ignore) {
  static double pik[K];
  int ans = 0;
  double pimax = log(0.0);
  for(int k=0; k<K; k++) {
    pik[k] = mcmclib_mvnorm_lpdf_compute_noinv(pi_hat[k], x);
    if(pik[k] > pimax) {
      pimax = pik[k];
      ans = k;
    }
  }
  return ans;
}

#include "ex4_target_distrib.c"

int main(int argc, char** argv) {
  /*set starting guess metropolis covariance matrices*/
  gsl_matrix* Sigma_local[K];
  for(int k=0; k < K; k++){
    Sigma_local[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_local[k]);
    gsl_matrix_scale(Sigma_local[k], SIGMA0_LOCAL * SCALING_FACTOR * V0[k]);
  }
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, SIGMA0_GLOBAL * SCALING_FACTOR);

  /**********************/
  /*INIT data structures*/
  /**********************/
  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*current chain value*/
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);
  /*init target distribution data*/
  target_distrib_init();
  /*init EM algorithm data*/
  gsl_matrix* P_hat = gsl_matrix_alloc(N, K);
  gsl_vector* w_hat = gsl_vector_alloc(K);
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], gsl_rng_uniform(r) * MU0);
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
    pi_hat[k] = mcmclib_mvnorm_lpdf_alloc(mu_hat[k], Sigma_hat[k]->data);
  }

  /*alloc a new RAPT sampler*/
  mcmclib_rapt* sampler =
    mcmclib_rapt_alloc(r, /*RNG*/
		       target_logdensity, NULL,
		       x, T0, /*current chain value, burnin*/
		       Sigma_zero, K, Sigma_local, /*startng guesses*/
		       which_region, NULL);
  /*set custom proposal weights*/
  for(int k=0; k < K; k++) {
    for(int k1=0; k1 < K; k1++)
      gsl_matrix_set(sampler->lambda, k, k1, k==k1 ? 0.5 : 0.0);
    gsl_matrix_set(sampler->lambda, k, K, 0.5);
  }

  /*main MCMC loop*/
  int naccept=0;
  gsl_matrix* naccept_m = gsl_matrix_alloc(K, K+1);
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  for(int n=0; n<N; n++) {
    mcmclib_rapt_update(sampler);
    mcmclib_rapt_update_proposals(sampler);
    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x);
    if(((n+1) % N0)==0) {
      gsl_matrix_view Xv = gsl_matrix_submatrix(X, 0, 0, n+1, DIM);
      mcmclib_mixem_fit(&(Xv.matrix), K, mu_hat, Sigma_hat, P_hat, w_hat, 1);
      for(int k=0; k<K; k++)
	mcmclib_mvnorm_lpdf_inverse(pi_hat[k]);
    }

    if(sampler->accepted) {
      naccept++;
      gsl_matrix_set(naccept_m, sampler->which_region_x, sampler->which_proposal,
		     gsl_matrix_get(naccept_m, sampler->which_region_x, sampler->which_proposal) + 1.0);
    }
  }
  /*print out final summary*/
  printf("Acceptance rate: %f\n",
	 (gsl_matrix_get(naccept_m, 0, 0) + gsl_matrix_get(naccept_m, 1, 1)) /
	 (gsl_matrix_get(sampler->visits, 0, 0) + gsl_matrix_get(sampler->visits, 1, 1)));

  /*store sampled values*/
  FILE* out_X = fopen("ex4_X.csv", "w");
  gsl_matrix_fprintf(out_X, X, "%f");
  fclose(out_X);
  /*store means and variances estimates*/
  FILE* out_mu = fopen("ex4_mu_hat.csv", "w");
  FILE* out_Sigma = fopen("ex4_Sigma_hat.csv", "w");
  for(int k=0; k<K; k++){
    gsl_vector_fprintf(out_mu, mu_hat[k], "%f");
    gsl_matrix_fprintf(out_Sigma, Sigma_hat[k], "%f");
  }
  fclose(out_mu);
  fclose(out_Sigma);
  /*store mixture weights estimates*/
  FILE* out_w = fopen("ex4_w_hat.csv", "w");
  gsl_vector_fprintf(out_w, w_hat, "%f");
  fclose(out_w);

  gsl_matrix_free(naccept_m);
  target_distrib_free();
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  for(int k=0; k<K; k++)
    gsl_matrix_free(Sigma_local[k]);
  gsl_matrix_free(Sigma_zero);
}
