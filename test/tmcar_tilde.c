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
  fflush(stdout);
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
  /************************/
  /*check Givens rotations*/
  /************************/
  gsl_vector* alpha = gsl_vector_alloc(3);
  gsl_vector_set(alpha, 0, -0.5);
  gsl_vector_set(alpha, 1, 0.1);
  gsl_vector_set(alpha, 2, 0.5);
  gsl_matrix* A = gsl_matrix_alloc(3, 3);
  mcmclib_Givens_rotations(A, alpha);
  gsl_matrix* A1 = gsl_matrix_alloc(3, 3);
  gsl_matrix_memcpy(A1, A);
  mcmclib_matrix_inverse(A1);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      assert(check_dequal(gsl_matrix_get(A, i, j), gsl_matrix_get(A1, j, i)));
  gsl_matrix_free(A1);
  gsl_matrix_free(A);
  gsl_vector_free(alpha);

  /************************/
  /*check MCAR distrib.   */
  /************************/
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

  for(int i=0; i<(P-1); i++)
    for(int j=i+1; j<P; j++) {
      gsl_matrix_set(p->Gamma, i, j, 0.6);
      gsl_matrix_set(p->Gamma, j, i, 0.6);
    }

  mcmclib_mcar_tilde_lpdf_update_blocks(p);
  mcmclib_mcar_tilde_lpdf_update_vcov(p);
  assert(check_dequal(gsl_matrix_get(p->vcov, 0, 0), 1.363636));
  assert(check_dequal(gsl_matrix_get(p->vcov, 1, 3), 0.290909));
  assert(check_dequal(gsl_matrix_get(p->vcov, 4, 1), 0.109091));
  assert(check_dequal(gsl_matrix_get(p->vcov, 5, 3), 0.145455));

  gsl_vector* x = gsl_vector_alloc(N*P);
  gsl_vector_set_all(x, 1.0);
  mcmclib_mcar_tilde_lpdf_compute(p, x);
  gsl_vector_free(x);

  mcmclib_mcar_tilde_lpdf_free(p);
  gsl_matrix_free(W);
  gsl_vector_free(mu);
}
