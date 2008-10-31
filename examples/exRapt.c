/**RAPT example*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>

#define OUTPUT_FILE "data.csv"
/*chain length*/
#define N 100000
/*burn in length*/
#define T0 200
/*initial covariance guess*/
#define V0 1.0

/*state space dimension*/
#define DIM 3

/*target distribution: uniform in the unit cube*/
double target_logdensity(void* ignore, gsl_vector* x) {
  int d = x->size;
  for(int i=0; i<d; i++) {
    double xi = gsl_vector_get(x, i);
    if((xi < 0.0) || (xi > 1.0))
      return log(0.0);
  }
  return 0;
}

/*boundary function*/
int which_region(gsl_vector* x, void* ignore) {
  if(gsl_vector_get(x, 0) < 0.5)
    return 0;
  else
    return 1;
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
  gsl_vector_set_all(x, 0.5);

  /*alloc a new RAPT sampler*/
  mcmclib_rapt* sampler = mcmclib_rapt_alloc(r,
					     target_logdensity, NULL,
					     x, T0, Sigma_zero,
					     2, Sigma_local,
					     which_region, NULL);

  /*open output file*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  /*print out csv header*/
  fprintf(out, "x1, x2, x3\n");
  
  /*main MCMC loop*/
  for(int i=0; i<N; i++) {
    mcmclib_rapt_update(sampler);
    for(int j=0; j<(d-1); j++)
      fprintf(out, "%f, ", gsl_vector_get(x, j));
    fprintf(out, "%f\n", gsl_vector_get(x, d-1));
  }
  
  fclose(out);
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  gsl_matrix_free(Sigma_zero);
}
