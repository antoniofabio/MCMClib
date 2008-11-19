/**
INCA example
 Takes 5.4sec on 800Mhz AMD Linux
*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_inca.h>

/*target space dimension*/
#define DIM 3
/*total number of iterations*/
#define N 100000
/*HCL burn in*/
#define T0 20
/*starting variance guess*/
#define V0 0.1
/*number of parallel chains to run*/
#define K 3

/*target distribution: uniform in the unit cube (side length=10)*/
double target_logdensity(void* ignore, gsl_vector* x) {
  int d = x->size;
  for(int i=0; i<d; i++) {
    double xi = gsl_vector_get(x, i);
    if((xi < -5.0) || (xi > 5.0))
      return log(0.0);
  }
  return 0;
}

int main(int argc, char** argv) {
  int d = DIM;
  /*set starting guess covariance matrix*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, V0);

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*lets go sample from K different chains in parallel*/
  mcmclib_gauss_inca_pool* pool = mcmclib_gauss_inca_pool_alloc(Sigma_zero, T0, K);
  gsl_vector* xx[K];
  mcmclib_gauss_inca* sampler[K];
  for(int k=0; k<K; k++) {
    sampler[k] = mcmclib_gauss_inca_alloc(pool);
    xx[k] = gsl_vector_alloc(d);
    /*set starting value at random*/
    for(int j=0; j<d; j++)
      gsl_vector_set(xx[k], j, (gsl_rng_uniform(r) * 10.0) - 5.0);
  }

  /*open output files*/
  FILE* out[K];
  char filename[512];
  for(int i=0; i<K; i++) {
    sprintf(filename, "data_%d.csv", i);
    out[i] = fopen(filename, "w");
    /*print out csv header*/
    fprintf(out[i], "x1, x2, x3\n");
  }

  /*main MCMC loop: for each time step iterate troughout the K chains*/
  for(int i=0; i<N; i++) for(int k=0; k<K; k++) {
      mcmclib_gauss_inca_update(sampler[k], r, target_logdensity, xx[k], NULL);
      for(int j=0; j<(d-1); j++)
	fprintf(out[k], "%f, ", gsl_vector_get(xx[k], j));
      fprintf(out[k], "%f\n", gsl_vector_get(xx[k], d-1));
    }

  /*release system resources*/
  for(int k=0; k<K; k++) {
    fclose(out[k]);
    mcmclib_gauss_inca_free(sampler[k]);
    gsl_vector_free(xx[k]);
  }
  mcmclib_gauss_inca_pool_free(pool);
  gsl_rng_free(r);
  gsl_matrix_free(Sigma_zero);
  return 0;
}
