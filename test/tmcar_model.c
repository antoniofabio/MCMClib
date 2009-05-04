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

void matrix_printf(gsl_matrix* m) {
  int p = m->size1;
  int q = m->size2;
  for(int i=0; i<p; i++) {
    for(int j=0; j<q; j++)
      printf("\t%f", gsl_matrix_get(m, i, j));
    printf("\n");
  }
  printf("\n");
  fflush(stdout);
}

#define P 3
#define N 4

mcmclib_mcar_model* p;

/* declare as near regions 'i' and 'j' */
#define DECL_AD(i, j) if(1) {			\
    gsl_matrix_set(W, i, j, 1.0);		\
    gsl_matrix_set(W, j, i, 1.0);		\
  }

double lpdf_alpha1(double s) {
  gsl_vector* alpha_h = gsl_vector_alloc(P * (P-1) / 2);
  gsl_vector_set_all(alpha_h, s);
  double ans = mcmclib_mcar_model_alpha1_lpdf(p, alpha_h);
  printf("%f -> %f\n", s, ans);
  gsl_vector_free(alpha_h);
  return ans;
}

double lpdf_alpha2(double s) {
  gsl_vector* alpha_h = gsl_vector_alloc(P * (P-1) / 2);
  gsl_vector_set_all(alpha_h, s);
  double ans = mcmclib_mcar_model_alpha2_lpdf(p, alpha_h);
  printf("%f -> %f\n", s, ans);
  gsl_vector_free(alpha_h);
  return ans;
}

double lpdf_sigma(double s) {
  gsl_vector* sigma = gsl_vector_alloc(P);
  gsl_vector_set_all(sigma, s);
  gsl_vector_set_all(sigma, 0.0);
  gsl_vector_set(sigma, 0, s);
  double ans = mcmclib_mcar_model_sigma_lpdf(p, sigma);
  printf("%f -> %f\n", s, ans);
  gsl_vector_free(sigma);
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

  lpdf_sigma(0.68);
  exit(1);

  gsl_vector* alpha_h = gsl_vector_alloc(P * (P-1) / 2);
  gsl_vector_set_all(alpha_h, 10.0);
  printf("alpha1:\n");
  double l1 = lpdf_alpha1(10.0);
  lpdf_alpha1(5.0);
  lpdf_alpha1(-10.0);
  assert(l1 == lpdf_alpha1(10.0));

  printf("alpha2:\n");
  l1 = lpdf_alpha2(10.0);
  lpdf_alpha2(5.0);
  lpdf_alpha2(-10.0);
  assert(l1 == lpdf_alpha2(10.0));

  printf("sigma:\n");
  for(double s=0.625; s<0.673; s+=0.005)
    lpdf_sigma(s);

  gsl_vector* sigma = gsl_vector_alloc(P);
  gsl_vector_set_all(sigma, 0.0);
  gsl_vector_set(sigma, 1, -1.0);
  l1 = mcmclib_mcar_model_sigma_lpdf(p, sigma);
  gsl_vector_set_all(sigma, 0.0);
  gsl_vector_set(sigma, 0, -1.0);
  double l2 = mcmclib_mcar_model_sigma_lpdf(p, sigma);
  printf("%f == %f ?\n", l1, l2);
  assert(gsl_finite(l1));
  assert(gsl_finite(l2));
  /* distrib. should not change for permutations of sigma: */
  assert(l1 == l2);
  gsl_vector_free(sigma);

  gsl_vector* alphasigma = gsl_vector_alloc(P * (P-1) / 2 + P);
  gsl_vector_set_zero(alphasigma);
  printf("%f -> %f\n", 0.0, mcmclib_mcar_model_alphasigma_lpdf(p, alphasigma));
  gsl_vector_free(alphasigma);

  mcmclib_mcar_model_free(p);
  gsl_vector_free(e);
  mcmclib_mcar_tilde_lpdf_free(llik);
  gsl_matrix_free(W);
}
