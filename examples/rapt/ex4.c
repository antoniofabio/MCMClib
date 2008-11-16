/**RAPT example 4: mixture-based adaptive boundary*/
#include <stdio.h>
#include <memory.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_mrw.h>
#include <rapt.h>
#include <mixem_rec.h>
#include <mixem.h>
#include <mvnorm.h>

/** Following arguments can be (almost) freely customized */
/*chain length*/
static int N = 100000;
/*state space dimension*/
static int DIM = 2;
/*absolute local mean value*/
static double MU0 = 1.5;
/*local variances multipliers*/
static double V0[] = {1.0, 1.0};
/*local pairwise correlations*/
static double RHO[] = {0.0, 0.0};
/*mixture proportion of component 1*/
static double BETA = 0.5;
/*burn in length*/
#define T0 ((DIM + DIM * (DIM-1) / 2) * 200)
/*update boundary every N0 iterations*/
#define N0 (N / 100)
/*scaling factor*/
#define SCALING_FACTOR (2.38 * 2.38 / (double) DIM)
/*starting local variance guess as scaling factor w.r.t. true value*/
#define SIGMA0_LOCAL 0.25
/*starting global variance guess*/
#define SIGMA0_GLOBAL pow(MU0 * 2.0, 2.0)
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
    pik[k] = mcmclib_mvnorm_lpdf_compute(pi_hat[k], x);
    if(pik[k] > pimax) {
      pimax = pik[k];
      ans = k;
    }
  }
  return ans;
}

#include "ex4_target_distrib.c"

int main(int argc, char** argv) {
  /*******************************/
  /*read input data from cmd line*/
  /*******************************/
  if(argc > 9) {
    printf("Passed %d arguments. Usage:\n", argc-1);
    printf("%s N DIM S1 S2 RHO1 RHO2 BETA\n", argv[0]);
    exit(1);
  }
#define scanfif(i, fmt, target) if(argc > (i)) sscanf(argv[(i)], fmt, target)
  scanfif(1, "%d", &N);
  scanfif(2, "%d", &DIM);
  scanfif(3, "%lf", V0);
  scanfif(4, "%lf", V0+1);
  scanfif(5, "%lf", RHO);
  scanfif(6, "%lf", RHO+1);
  scanfif(7, "%lf", &BETA);
  scanfif(8, "%lf", &MU0);

  printf("=Simulation settings=\n");
  printf("N = %d\tDIM = %d\n", N, DIM);
  printf("S1 = %f\tS2 = %f\n", V0[0], V0[1]);
  printf("RHO1 = %f\tRHO2 = %f\n", RHO[0], RHO[1]);
  printf("BETA = %f, MU0 = %f\n", BETA, MU0);
  printf("T0 = %d\tN0 = %d\n", T0, N0);
  printf("=====================\n");
  /*******************************/

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

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*current chain value*/
  gsl_vector* x = gsl_vector_alloc(DIM);
  gsl_vector_set_all(x, 0.0);
  /*init target distribution data*/
  target_distrib_init();
  /*init EM algorithm data*/
  gsl_vector* beta_hat = gsl_vector_alloc(K);
  for(int k=0; k<K; k++) {
    mu_hat[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu_hat[k], gsl_rng_uniform(r) * MU0 * pow(-1.0, k+1));
    Sigma_hat[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma_hat[k]);
    gsl_matrix_scale(Sigma_hat[k], SIGMA0_LOCAL);
    pi_hat[k] = mcmclib_mvnorm_lpdf_alloc(mu_hat[k], Sigma_hat[k]->data);
  }

  /*alloc a new RAPT sampler*/
  mcmclib_rapt* sampler =
    mcmclib_rapt_alloc(r, /*RNG*/
		       target_logdensity, NULL, /*target distribution*/
		       x, T0, /*current chain value, burnin*/
		       Sigma_zero, K, Sigma_local, /*starting guesses*/
		       which_region, NULL); /*boundary function*/
  /*set custom proposal weights*/
  for(int k=0; k < K; k++) {
    for(int k1=0; k1 < K; k1++)
      gsl_matrix_set(sampler->lambda, k, k1, k==k1 ? 0.5 : 0.0);
    gsl_matrix_set(sampler->lambda, k, K, 0.5);
  }

  /*main MCMC loop*/
  int naccept=0; /*number of acceptances*/
  gsl_matrix* naccept_m = gsl_matrix_alloc(K, K+1); /*# accept. x region & proposal*/
  gsl_matrix* X = gsl_matrix_alloc(N, DIM); /*matrix of all sampled values*/
  mcmclib_mixem_rec* m = mcmclib_mixem_rec_alloc(mu_hat, Sigma_hat, beta_hat);
  for(int n=0; n<N; n++) {
    /*update chain value*/
    mcmclib_rapt_update(sampler);
    /*store new value in matrix X*/
    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x);

    /*accumulate data in mixture fitter object*/
    mcmclib_mixem_rec_add(m, x);
    /*update boundary estimation and rapt proposal variances*/
    if((n>T0) && (((n+1) % N0)==0)) {
      mcmclib_rapt_update_proposals(sampler);
      mcmclib_mixem_rec_update(m);
    }

    /*update acceptance rate informations*/
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
  double td = fabs(gsl_vector_get(beta_hat, 0) - BETA);
  printf("Distance between true and estimated boundary: %f\n", td);

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
  FILE* out_beta = fopen("ex4_beta_hat.csv", "w");
  gsl_vector_fprintf(out_beta, beta_hat, "%f");
  fclose(out_beta);

  gsl_matrix_free(naccept_m);
  target_distrib_free();
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  for(int k=0; k<K; k++)
    gsl_matrix_free(Sigma_local[k]);
  gsl_matrix_free(Sigma_zero);
}
