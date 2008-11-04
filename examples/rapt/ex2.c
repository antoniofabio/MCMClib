/**RAPT example 2: moving the boundary*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>

#include "ex2_extra.c"

#define OUTPUT_FILE "ex2_out.csv"
#define EXTRA_OUTPUT_FILE "ex2_extra_out.csv"
/*chain blocks length*/
#define B0 1e5
/*chain length as number of blocks*/
#define N 10
/*burn in length as number of its.*/
#define T0 200
/*initial variance guess*/
#define V0 0.5

/*state space dimension*/
#define DIM 2

/*target distribution: uniform with 'V' shape*/
double target_logdensity(void* ignore, gsl_vector* vx) {
  double x = gsl_vector_get(vx, 0);
  double y = gsl_vector_get(vx, 1);
  double width = 0.1;
  if(fabs(x) > 1.0)
    return log(0.0);
  if(y < -(width/2.0) )
    return log(0.0);
  if(fabs(y - fabs(x)) > width)
    return log(0.0);
  return 0.0;
}

/*boundary function*/
int which_region(gsl_vector* x, void* in_th) {
  double th = *((double*) in_th);
  if(gsl_vector_get(x, 0) < th)
    return 0;
  else
    return 1;
}

static void print_vector(FILE* stream, gsl_vector* x) {
  for(int i=0; i< (x->size - 1); i++)
    fprintf(stream, "%.3f, ", gsl_vector_get(x, i));
  fprintf(stream, "%.3f", gsl_vector_get(x, x->size -1));
}

int main(int argc, char** argv) {
  int d = DIM;
  /*set starting guess covariance matrix*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, V0);
  gsl_matrix** Sigma_local = (gsl_matrix**) malloc(2 * sizeof(gsl_matrix*));
  for(int k=0; k<2; k++){
    Sigma_local[k] = gsl_matrix_alloc(d,d);
    gsl_matrix_memcpy(Sigma_local[k], Sigma_zero);
  }

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*current chain value*/
  gsl_vector* x = gsl_vector_alloc(d);
  gsl_vector_set_all(x, 0.0);
  /*current threshold value*/
  double th = -0.5;

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

  /*open output file*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  /*print out csv header*/
  for(int j=0; j<d; j++)
    fprintf(out, "x%d, ", j);
  fprintf(out, "proposal\n");

  /*alloc block info data*/
  block_info* bi = block_info_alloc(sampler, B0);

  /*main MCMC loop*/
  double score=1e6;
  double newscore;
  double oldth=th;
  for(int b=0; b<N; b++) {
    for(int nb=0; nb< B0; nb++) {
      mcmclib_rapt_update(sampler);
      mcmclib_rapt_update_proposals(sampler);
      /*      print_vector(out, x);
	      fprintf(out, ", %d\n", sampler->which_proposal);*/
      block_info_update(bi);
    }
    newscore = block_info_score(bi);
    printf("th = %f; ", th);
    print_vector(stdout, bi->den);
    printf("; "); print_vector(stdout, bi->num);
    printf(" -> %f\n", score);
    score = newscore;
    oldth = th;
    //th += (0.1 / (b + 1.0)) * ( (gsl_rng_uniform(r) < 0.5) ? -1.0 : 1.0 );
    th += 0.1;
  }

  fclose(out);
  block_info_free(bi);
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  gsl_matrix_free(Sigma_zero);
}
