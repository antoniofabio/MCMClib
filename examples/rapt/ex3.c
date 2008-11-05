/**RAPT example 3: moving the boundary, high dimensionality*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>
#include <mvnorm.h>

#include "ex2_extra.c"

#define OUTPUT_FILE "ex3_out.csv"
#define EXTRA_OUTPUT_FILE "ex3_extra_out.csv"
/*chain blocks length*/
#define B0 ( 400 * (DIM + DIM * (DIM-1) / 2) )
/*chain length as number of blocks*/
#define N 100
/*burn in length as number of its.*/
#define T0 B0/2
/*initial variance guess*/
#define V0 1.0
/*initial threshold guess*/
#define TH0 0.0

/*state space dimension*/
#define DIM 100

/*definitions for the target distribution*/
#define ABS_MU 1.5
static gsl_vector* mu1;
static gsl_vector* mu2;
static gsl_matrix* Sigma1;
static gsl_matrix* Sigma2;
static mcmclib_mvnorm_lpdf* pi1;
static mcmclib_mvnorm_lpdf* pi2;
/*target log-distribution: mixture of two normals*/
double target_logdensity(void* ignore, gsl_vector* x) {
  return log(exp(mcmclib_mvnorm_lpdf_compute(pi1, x)) +
	     exp(mcmclib_mvnorm_lpdf_compute(pi2, x)));
}
static void target_distrib_init() {
  mu1 = gsl_vector_alloc(DIM);
  mu2 = gsl_vector_alloc(DIM);
  Sigma1 = gsl_matrix_alloc(DIM, DIM);
  Sigma2 = gsl_matrix_alloc(DIM, DIM);
  gsl_vector_set_all(mu1, ABS_MU);
  gsl_vector_set_all(mu2, - ABS_MU);

  gsl_matrix_set_identity(Sigma1);
  gsl_matrix_scale(Sigma1, V0);
  gsl_matrix_memcpy(Sigma2, Sigma1);
  /*  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_matrix* tmp1 = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix* tmp2 = gsl_matrix_alloc(DIM, DIM);
  for(int i=0; i<DIM; i++) for(int j=0; j<DIM; j++) {
      gsl_matrix_set(tmp1, i, j, gsl_ran_gaussian(r, 1.0));
      gsl_matrix_set(tmp2, i, j, gsl_ran_gaussian(r, 1.0));
    }
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp1, tmp1, 0.0, Sigma1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp2, tmp2, 0.0, Sigma2);
  gsl_matrix_free(tmp1);
  gsl_matrix_free(tmp2);
  gsl_rng_free(r);*/
  pi1 = mcmclib_mvnorm_lpdf_alloc(mu1, Sigma1->data);
  pi2 = mcmclib_mvnorm_lpdf_alloc(mu2, Sigma2->data);
}
static void target_distrib_free() {
  gsl_vector_free(mu1);
  gsl_vector_free(mu2);
  gsl_matrix_free(Sigma1);
  gsl_matrix_free(Sigma2);
  mcmclib_mvnorm_lpdf_free(pi1);
  mcmclib_mvnorm_lpdf_free(pi2);
}

/*boundary function*/
int which_region(gsl_vector* x, void* in_th) {
  double th = *((double*) in_th);
  double tmp = 0.0;
  for(int i=0; i<DIM; i++)
    tmp += gsl_vector_get(x, i);
  //  tmp = gsl_vector_get(x, 0);
  //  printf("sum(x) = %f\n", tmp);
  if(tmp < th)
    return 0;
  else
    return 1;
}

int main(int argc, char** argv) {
  int d = DIM;
  /*set starting guess covariance matrix*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, V0 * 2.38 * 2.38 / DIM);
  gsl_matrix** Sigma_local = (gsl_matrix**) malloc(2 * sizeof(gsl_matrix*));
  for(int k=0; k<2; k++){
    Sigma_local[k] = gsl_matrix_alloc(d,d);
    gsl_matrix_memcpy(Sigma_local[k], Sigma_zero);
  }
  gsl_matrix_scale(Sigma_zero, 1.0 * DIM);

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*current chain value*/
  gsl_vector* x = gsl_vector_alloc(d);
  gsl_vector_set_all(x, 0.0);
  /*starting threshold value*/
  double th = TH0;

  /*alloc a new RAPT sampler*/
  mcmclib_rapt* sampler = mcmclib_rapt_alloc(r,
					     target_logdensity, NULL,
					     x, T0, Sigma_zero,
					     2, Sigma_local,
					     which_region, &th);
  for(int k=0; k<2; k++) {
    for(int k1=0; k1<2; k1++)
      gsl_matrix_set(sampler->lambda, k, k1, k==k1 ? 0.5 : 0.0);
    gsl_matrix_set(sampler->lambda, k, 2, 0.5);
  }

  /*open output files*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  FILE* out_extra = fopen(EXTRA_OUTPUT_FILE, "w");
  fprintf(out_extra, "th, score\n");

  /*alloc block info data*/
  block_info* bi = block_info_alloc(sampler, B0);

  target_distrib_init();

  /*main MCMC loop*/
  double score=1e6;
  double newscore;
  double oldth=th;
  for(int b=0; b<N; b++) {
    for(int nb=0; nb< B0; nb++) {
      mcmclib_rapt_update(sampler);
      mcmclib_rapt_update_proposals(sampler);
      block_info_update(bi);
      for(int i=0; i<DIM; i++)
	fprintf(out, "%f ", gsl_vector_get(x, i));
      fprintf(out, "\n");
    }
    newscore = block_info_score(bi);
    //printf("th = %f -> %f\n", th, newscore);
    fprintf(out_extra, "%f, %f\n", th, newscore);
    if(newscore > score) {
      //printf("oops! go back to %f...\n", oldth);
      th = oldth;
      newscore = score;
    }
    /*    printf(".");
    if((b % 10) == 0)
      printf("\n%d\n", b);
      fflush(stdout);*/
    score = newscore;
    oldth = th;
    /*    double width = 0.2 / pow((10 * (b + 1)/ (double) N), 2);
    th += gsl_rng_uniform(r) * width - width/2.0;
    if(th < -ABS_MU) th = -ABS_MU;
    else if(th > ABS_MU) th = ABS_MU;*/
  }

  target_distrib_free();
  fclose(out_extra);
  block_info_free(bi);
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  gsl_matrix_free(Sigma_zero);
}
