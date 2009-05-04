/**Test spatial normal distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mcar_tilde.h>

#define TOL 1e-6
int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

void mprint(gsl_matrix* A) {
  int n = A->size1;
  int p = A->size2;
  for(int i=0; i<n; i++) {
    for(int j=0; j<p; j++)
      printf("%.3f, ", gsl_matrix_get(A, i, j));
    printf("\n");
  }
}

int is_pos_def(gsl_matrix* A) {
  if(A->size1 != A->size2)
    return 0;
  int n = A->size1;
  for(int i=0; i<(n-1); i++) {
    if(gsl_matrix_get(A, i, i) <= 0.0)
      return 0;
    for(int j=(i+1); j<n; j++)
      if(fabs(gsl_matrix_get(A, i, j) - gsl_matrix_get(A, j, i)) >= TOL)
	return 0;
  }
  return 1;
}

#define P 2
#define N 3

mcmclib_mcar_tilde_lpdf* p;

/* declare as near regions 'i' and 'j' */
#define DECL_AD(i, j) if(1) {			\
    gsl_matrix_set(W, i, j, 1.0);		\
    gsl_matrix_set(W, j, i, 1.0);		\
  }

gsl_vector* x;
static double lpdf(double s) {
  gsl_vector_set_all(x, s);
  return mcmclib_mcar_tilde_lpdf_compute(p, x);
}

int main(int argc, char** argv) {
  gsl_vector* mu = gsl_vector_alloc(N * P);
  gsl_vector_set_zero(mu);
  gsl_matrix* W = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(W);
  for(int i=0; i<(N-1); i++)
    DECL_AD(i, i+1);
  
  p = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  for(int i=0; i<N; i++) {
    int count = 0;
    for(int j=0; j<N; j++)
      count += gsl_matrix_get(W, i, j) == 1.0;
    assert(gsl_vector_get(p->m, i) == (double) count);
  }

  gsl_vector_set_all(p->alphasigmag, 0.0);
  gsl_vector_set_all(p->alpha12sigma, -1.0);

  mcmclib_mcar_tilde_lpdf_update_blocks(p);
  printf("Gamma:\n");
  mprint(p->Gamma);
  printf("B_tilde:\n");
  mprint(p->B_tilde);
  printf("Blocks:\n");
  mprint(p->vcov);
  mcmclib_mcar_tilde_lpdf_update_vcov(p);
  printf("VCOV:\n");
  mprint(p->vcov);
  assert(is_pos_def(p->vcov));

  gsl_vector_set_all(p->alpha12sigma, 0.7);
  mcmclib_mcar_tilde_lpdf_update_blocks(p);
  printf("B_tilde(0.7):\n");
  mprint(p->B_tilde);
  printf("Blocks(0.7):\n");
  mprint(p->vcov);
  mcmclib_mcar_tilde_lpdf_update_blocks(p);

  mcmclib_mcar_tilde_lpdf_update_vcov(p);
  printf("VCOV(0.7):\n");
  mprint(p->vcov);
  assert(is_pos_def(p->vcov));

  x = gsl_vector_alloc(N*P);
  printf("%f -> %f\n", 0.0, lpdf(0.0));
  printf("%f -> %f\n", 1.0, lpdf(1.0));
  gsl_vector_set_all(p->alpha12sigma, 0.7);
  printf("%f -> %f\n", 1.0, lpdf(1.0));

  gsl_vector_free(x);

  mcmclib_mcar_tilde_lpdf_free(p);
  gsl_matrix_free(W);
  gsl_vector_free(mu);
}
