/**RAPT example 2*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>

/*trace program execution to stdout?*/
//#define TRACE_ME

#define OUTPUT_FILE "data_rapt2.csv"
/*chain length*/
#define N 1e5
/*burn in length*/
#define T0 200
/*initial covariance guess*/
#define V0 0.5

/*state space dimension*/
#define DIM 2

/*target distribution: 'V' shape*/
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
int which_region(gsl_vector* x, void* ignore) {
  if(gsl_vector_get(x, 0) < -0.5)
    return 0;
  else
    return 1;
}

static void print_vector(FILE* stream, gsl_vector* x) {
  for(int i=0; i< (x->size - 1); i++)
    fprintf(stream, "%.3f, ", gsl_vector_get(x, i));
  fprintf(stream, "%.3f", gsl_vector_get(x, x->size -1));
}

static void dump_vector(gsl_vector* x) {
  printf("(");
  print_vector(stdout, x);
  printf(")\n");
}
static void dump_matrix(gsl_matrix* x) {
  for(int j=0; j< x->size1; j++) {
    gsl_vector_view rv = gsl_matrix_row(x, j);
    dump_vector(&(rv.vector));
  }
}

/*dump current sampler state to stdout*/
static void dump_rapt(mcmclib_rapt* s) {
  printf("#current_x: ");  dump_vector(s->current_x);
  printf("#old: "); dump_vector(s->old);
  printf("#accepted: %d\n", s->accepted);
  printf("#which_proposal: %d\n", s->which_proposal);
  printf("#ntries: "), dump_vector(s->ntries);
  printf("#sigma_whole:\n"); dump_matrix(s->sigma_whole);
  printf("#sigma_local[...]:\n");
  for(int k=0; k< s->K; k++) {
    printf(" [%d]:\n", k);
    dump_matrix(s->sigma_local[k]);
  }
  printf("#t: %d\n", s->t);
  printf("#global_mean: "); dump_vector(s->global_mean);
  printf("#global_variance:\n"); dump_matrix(s->global_variance);
  printf("#means[...]:\n");
  for(int k=0; k< s->K; k++) {
    printf(" [%d]: ", k);
    dump_vector(s->means[k]);
  }  
  printf("#variances[...]:\n");
  for(int k=0; k< s->K; k++) {
    printf(" [%d]:\n", k);
    dump_matrix(s->variances[k]);
  }
  printf("#n: "); dump_vector(s->n);
  printf("#visits:\n"); dump_matrix(s->visits);
  printf("#jd:\n"); dump_matrix(s->jd);
  printf("#lambda:\n"); dump_matrix(s->lambda);
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

  /*alloc a new RAPT sampler*/
  mcmclib_rapt* sampler = mcmclib_rapt_alloc(r,
					     target_logdensity, NULL,
					     x, T0, Sigma_zero,
					     2, Sigma_local,
					     which_region, NULL);

  /*open output file*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  /*print out csv header*/
  for(int j=0; j<d; j++)
    fprintf(out, "x%d, ", j);
  fprintf(out, "proposal\n");

  /*open extra output file*/
  FILE* out_extra = fopen("out_rapt_extra.csv", "w");
  /*print out csv header*/
  for(int j=0; j<d; j++)
    fprintf(out_extra, "x%d, ", j);
  fprintf(out_extra, "proposal, ntries0, ntries1, ntries2\n");

  /*main MCMC loop*/
  for(int i=0; i<N; i++) {
#ifdef TRACE_ME
    printf("\n-----\niteration %d\nsampler:\n", i);
    dump_rapt(sampler);
#endif
    mcmclib_rapt_update(sampler);
    mcmclib_rapt_update_proposals(sampler);
    print_vector(out, x);
    fprintf(out, "%d\n", sampler->which_proposal);
    if(sampler->accepted) {
      print_vector(out_extra, sampler->old);
      fprintf(out_extra, ", %d, ", sampler->which_proposal);
      print_vector(out_extra, sampler->ntries);
      fprintf(out_extra, "\n");
    }
  }
  
  fclose(out_extra);
  fclose(out);
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  gsl_matrix_free(Sigma_zero);
}
