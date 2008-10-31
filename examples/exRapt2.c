/**RAPT example 2*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <rapt.h>

#define OUTPUT_FILE "data_rapt2.csv"
/*chain length*/
#define N 100000
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
  if(gsl_vector_get(x, 0) < 0.0)
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
  for(int j=0; j<(d-1); j++)
    fprintf(out, "x%d, ", j);
  fprintf(out, "x%d\n", d-1);
  
  /*main MCMC loop*/
  for(int i=0; i<N; i++) {
    mcmclib_rapt_update(sampler);
    /*    if((i % 10) == 0)
	  mcmclib_rapt_update_lambda(sampler);*/
    for(int j=0; j<(d-1); j++)
      fprintf(out, "%f, ", gsl_vector_get(x, j));
    fprintf(out, "%f\n", gsl_vector_get(x, d-1));
  }
  
  fclose(out);
  gsl_rng_free(r);
  mcmclib_rapt_free(sampler);
  gsl_matrix_free(Sigma_zero);
}
