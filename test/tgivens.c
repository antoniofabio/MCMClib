/** Test Givens rotations */
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mcar_tilde.h>
#include <givens.h>

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

#define P 2
#define N 3

int main(int argc, char** argv) {
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
}
