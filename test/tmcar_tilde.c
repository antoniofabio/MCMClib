/**Test spatial normal distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mcar_tilde.h>

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

void matrix_printf(gsl_matrix* m) {
  int p = m->size1;
  int q = m->size2;
  for(int i=0; i<p; i++) {
    for(int j=0; j<q; j++)
      printf("\t%f", gsl_matrix_get(m, i, j));
    printf("\n");
  }
  printf("\n");
}

#define P 2
#define N 3

mcmclib_mcar_tilde_lpdf* p;

/* declare as near regions 'i' and 'j' */
#define DECL_AD(i, j) if(1) {			\
    gsl_matrix_set(W, i, j, 1.0);		\
    gsl_matrix_set(W, j, i, 1.0);		\
  }

int main(int argc, char** argv) {
  gsl_vector* mu = gsl_vector_alloc(N * P);
  gsl_vector_set_zero(mu);
  gsl_matrix* W = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(W);
  for(int i=0; i<(N-1); i++)
    DECL_AD(i, i+1);
  
  p = mcmclib_mcar_tilde_lpdf_alloc(P, N, W);
  mcmclib_mcar_tilde_lpdf_update_vcov(p);
  matrix_printf(p->vcov);

  mcmclib_mcar_tilde_lpdf_free(p);
  gsl_matrix_free(W);
}
