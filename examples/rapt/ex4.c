/**RAPT example 4: mixture-based adaptive boundary*/
#include <stdio.h>
#include <memory.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>
#include <mixem.h>
#include <mvnorm.h>

#define OUTPUT_FILE "ex3_out.csv"
#define EXTRA_OUTPUT_FILE "ex3_extra_out.csv"

/*burn in length*/
#define T0 ((DIM + DIM * (DIM-1) / 2) * 100)
/*initial variance guess*/
#define V0 1.0
/*update boundary every N0 iterations*/
static int N0 = 10000;


/*input arguments*/
static int N = 100000;
static int DIM = 2;
static double MU0 = 1.5;
static void scan_arguments(int argc, char** argv) {
  if(argc==1) return;
  sscanf(argv[1], "%d", &N);
  if(argc==2) return;
  sscanf(argv[2], "%d", &DIM);
  if(argc == 3) return;
  sscanf(argv[3], "%lf", &MU0);
}

#define K 2
static gsl_vector* mu_hat[K];
static gsl_matrix* Sigma_hat[K];
static mcmclib_mvnorm_lpdf* pi_hat[K];
/*boundary function*/
int which_region(gsl_vector* x, void* ignore) {
  static double pik[K];
  int ans = 0;
  double pimax = log(0.0);
  for(int k=0; k<K; k++) {
    pik[k] = mcmclib_mvnorm_lpdf_compute_nochol(pi_hat[k], x);
    if(pik[k] > pimax) {
      pimax = pik[k];
      ans = k;
    }
  }
  return ans;
}

#include "ex3_target_distrib.c"

int main(int argc, char** argv) {
  /*scan input arguments*/
  if(argc > 4) {
    fprintf(stderr, "expecting arguments: N(=%d), dim(=%d), mu0(=%f)\n", N, DIM, MU0);
    exit(1);
  }
  scan_arguments(argc, argv);

  /*set starting guess covariance matrices*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, V0 * 2.38 * 2.38 / DIM);
  gsl_matrix** Sigma_local = (gsl_matrix**) malloc(2 * sizeof(gsl_matrix*));
  for(int k=0; k<2; k++){
    Sigma_local[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_memcpy(Sigma_local[k], Sigma_zero);
    gsl_matrix_scale(Sigma_local[k], 0.25);
  }
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, pow(2.0 * MU0, 2.0) * 2.38 * 2.38 / DIM);

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
		       Sigma_zero, 2, Sigma_local, /*startng guesses*/
		       which_region, NULL);
  /*set custom proposal weights*/
  for(int k=0; k<2; k++) {
    for(int k1=0; k1<2; k1++)
      gsl_matrix_set(sampler->lambda, k, k1, k==k1 ? 0.5 : 0.0);
    gsl_matrix_set(sampler->lambda, k, 2, 0.5);
  }

  /*open output files, write headers*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  for(int i=0; i<DIM; i++)
    fprintf(out, "x%d, ", i);
  fprintf(out, "proposal\n");
  FILE* out_extra = fopen(EXTRA_OUTPUT_FILE, "w");
  for(int i=0; i<DIM; i++)
    fprintf(out_extra, "x%d, ", i);
  fprintf(out_extra, "ntries0, ntries1, ntries2, jump, proposal\n");

  /*main MCMC loop*/
  int naccept=0;
  gsl_matrix* naccept_m = gsl_matrix_alloc(K, K+1);
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  for(int n=0; n<N; n++) {

#ifndef NO_SAVE
    if(n>0 && sampler->accepted) {
      for(int i=0; i<DIM; i++)
	fprintf(out_extra, "%f, ", gsl_vector_get(sampler->old, i));
      fprintf(out_extra, "%f, %f, %f, %f, %d\n",
	      gsl_vector_get(sampler->ntries, 0),
	      gsl_vector_get(sampler->ntries, 1),
	      gsl_vector_get(sampler->ntries, 2),
	      sampler->last_jd, sampler->which_proposal);
    }
#endif

    mcmclib_rapt_update(sampler);
    mcmclib_rapt_update_proposals(sampler);
    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x);
    if(((n+1) % N0)==0) {
      gsl_matrix_view Xv = gsl_matrix_submatrix(X, 0, 0, n, DIM);
      mcmclib_mixem_fit(&(Xv.matrix), K, mu_hat, Sigma_hat, P_hat, w_hat, 10);
      for(int k=0; k<K; k++)
	mcmclib_mvnorm_lpdf_chol(pi_hat[k]);
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

  printf("Means estimates:\n");
  for(int k=0; k<K; k++){
    printf("mu[%d]:\n", k);
    gsl_vector_fprintf(stdout, mu_hat[k], "%f");
  }
  printf("Variances estimates:\n");
  for(int k=0; k<K; k++){
    printf("Sigma[%d]:\n", k);
    gsl_matrix_fprintf(stdout, Sigma_hat[k], "%f");
  }

  gsl_matrix_free(naccept_m);
  fclose(out_extra);
  fclose(out);
  target_distrib_free();
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  for(int k=0; k<2; k++)
    gsl_matrix_free(Sigma_local[k]);
  gsl_matrix_free(Sigma_zero);
}
