/**Test MCAR model */
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <mcar_model.h>

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

#define P 3
#define N 4

mcmclib_mcar_model* p;

/* declare as near regions 'i' and 'j' */
#define DECL_AD(i, j) if(1) {			\
    gsl_matrix_set(W, i, j, 1.0);		\
    gsl_matrix_set(W, j, i, 1.0);		\
  }

double lpdf_alpha12sigma(double s) {
  gsl_vector* alpha12sigma = gsl_vector_alloc(P * P);
  gsl_vector_set_all(alpha12sigma, s);
  double ans = mcmclib_mcar_model_alpha12sigma_lpdf(p, alpha12sigma);
  printf("%f -> %f\n", s, ans);
  gsl_vector_free(alpha12sigma);
  return ans;
}

int main(int argc, char** argv) {
  gsl_matrix* W = gsl_matrix_alloc(N, N);
  gsl_matrix_set_zero(W);
  for(int i=0; i<(N-1); i++)
    DECL_AD(i, i+1);

  mcmclib_mcar_tilde_lpdf* llik = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  gsl_vector* e = gsl_vector_alloc(N * P);
  gsl_vector_set_all(e, 2.0);
  p = mcmclib_mcar_model_alloc(llik, e);

  double l1 = lpdf_alpha12sigma(5.0);
  lpdf_alpha12sigma(-10.0);
  assert(l1 == lpdf_alpha12sigma(5.0));

  gsl_vector* alphasigma = gsl_vector_alloc(P * (P-1) / 2 + P);
  gsl_vector_set_zero(alphasigma);
  printf("%f -> %f\n", 0.0, mcmclib_mcar_model_alphasigma_lpdf(p, alphasigma));
  gsl_vector_free(alphasigma);

  mcmclib_mcar_model_free(p);
  gsl_vector_free(e);
  mcmclib_mcar_tilde_lpdf_free(llik);
  gsl_matrix_free(W);
}
