/**RAPT example 3: user-specified mu, boundary, dimensions*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>
#include <mvnorm.h>

#define OUTPUT_FILE "ex3_out.csv"
#define EXTRA_OUTPUT_FILE "ex3_extra_out.csv"

/*burn in length*/
#define T0 ((DIM + DIM * (DIM-1) / 2) * 100)
/*initial variance guess*/
#define V0 1.0

/*input arguments*/
static double TH;
static int N = 100000;
static int DIM = 2;
static double MU0 = 1.5;
static void scan_arguments(int argc, char** argv) {
  sscanf(argv[1], "%lf", &TH);
  if(argc==2) return;
  sscanf(argv[2], "%d", &N);
  if(argc==3) return;
  sscanf(argv[3], "%d", &DIM);
  if(argc == 4) return;
  sscanf(argv[4], "%lf", &MU0);
}

/*boundary function*/
int which_region(gsl_vector* x, void* ignore) {
  double tmp = 0.0;
  for(int i=0; i<DIM; i++)
    tmp += gsl_vector_get(x, i);
  if(tmp < TH)
    return 0;
  else
    return 1;
}

#include "ex3_target_distrib.c"

int main(int argc, char** argv) {
  /*scan input arguments*/
  if(argc < 2) {
    fprintf(stderr, "expecting arguments: th, N(=%d), dim(=%d), mu0(=%f)\n",
	    N, DIM, MU0);
    exit(1);
  }
  scan_arguments(argc, argv);
  int d = DIM;

  /*set starting guess covariance matrices*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, V0 * 2.38 * 2.38 / DIM);
  gsl_matrix** Sigma_local = (gsl_matrix**) malloc(2 * sizeof(gsl_matrix*));
  for(int k=0; k<2; k++){
    Sigma_local[k] = gsl_matrix_alloc(d,d);
    gsl_matrix_memcpy(Sigma_local[k], Sigma_zero);
    gsl_matrix_scale(Sigma_local[k], 0.25);
  }

  /**********************/
  /*INIT data structures*/
  /**********************/
  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*current chain value*/
  gsl_vector* x = gsl_vector_alloc(d);
  gsl_vector_set_all(x, 0.0);

  /*init target distribution data*/
  target_distrib_init();

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
  for(int n=0; n<N; n++) {
    if(n>0 && sampler->accepted) {
      for(int i=0; i<DIM; i++)
	fprintf(out_extra, "%f, ", gsl_vector_get(sampler->old, i));
      fprintf(out_extra, "%f, %f, %f, %f, %d\n",
	      gsl_vector_get(sampler->ntries, 0),
	      gsl_vector_get(sampler->ntries, 1),
	      gsl_vector_get(sampler->ntries, 2),
	      sampler->last_jd, sampler->which_proposal);
    }

    mcmclib_rapt_update(sampler);
    mcmclib_rapt_update_proposals(sampler);

    for(int i=0; i<DIM; i++)
      fprintf(out, "%f, ", gsl_vector_get(x, i));
    fprintf(out, "%d\n", sampler->which_proposal);
  }

  fclose(out_extra);
  fclose(out);
  target_distrib_free();
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  for(int k=0; k<2; k++)
    gsl_matrix_free(Sigma_local[k]);
  gsl_matrix_free(Sigma_zero);
}
