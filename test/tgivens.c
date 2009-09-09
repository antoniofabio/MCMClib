/** Test Givens rotations */
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <matrix.h>
#include <givens.h>

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static double backTransform(double a) {
  return log((M_PI_2 + a) / (M_PI_2 - a));
}

void backTransformVector(gsl_vector* v) {
  int n = v->size;
  for(int i=0; i<n; i++)
    gsl_vector_set(v, i, backTransform(gsl_vector_get(v, i)));
}

gsl_matrix* A;
gsl_vector* alphasigma;

static void singValues(gsl_matrix* A, gsl_vector* values) {
  int n = A->size1;
  gsl_vector* work = gsl_vector_alloc(n);
  gsl_matrix* A1 = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(A1, A);
  gsl_matrix* V = gsl_matrix_alloc(n, n);
  gsl_linalg_SV_decomp (A1, V, values, work);
  gsl_matrix_free(V);
  gsl_matrix_free(A1);
  gsl_vector_free(work);
}

static void sRepresentation(double s) {
  gsl_vector_set_all(alphasigma, s);
  mcmclib_Givens_representation(A, alphasigma);
  gsl_vector* vals = gsl_vector_alloc(3);
  singValues(A, vals);
  for(int i=0; i<3; i++)
    assert(check_dequal(gsl_vector_get(vals, i), exp(s)));
  gsl_linalg_cholesky_decomp(A);
  gsl_vector_free(vals);
}

static void sRepresentationAsymm(double s) {
  gsl_vector_set_all(alphasigma, s);
  gsl_vector_set(alphasigma, 0, s + 2.0);
  gsl_vector_set(alphasigma, 1, s + 2.0);
  gsl_vector_set(alphasigma, 2, s + 2.0);
  int offset = 6;
  gsl_vector_set(alphasigma, offset, s - 1.0);
  gsl_vector_set(alphasigma, offset+2, s + 1.0);
  mcmclib_Givens_representation_asymm(A, alphasigma);
  gsl_vector* values = gsl_vector_alloc(3);
  singValues(A, values);
  assert(check_dequal(gsl_vector_get(values, 0), exp(s + 1.0)));
  assert(check_dequal(gsl_vector_get(values, 1), exp(s)));
  assert(check_dequal(gsl_vector_get(values, 2), exp(s - 1.0)));
}

int main(int argc, char** argv) {
  gsl_vector* alpha = gsl_vector_alloc(3);
  gsl_vector_set(alpha, 0, -0.5);
  gsl_vector_set(alpha, 1, 0.1);
  gsl_vector_set(alpha, 2, 0.5);
  A = gsl_matrix_alloc(3, 3);
  mcmclib_Givens_rotations(A, alpha);
  gsl_matrix* A1 = gsl_matrix_alloc(3, 3);
  gsl_matrix_memcpy(A1, A);
  mcmclib_matrix_inverse(A1);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      assert(check_dequal(gsl_matrix_get(A, i, j), gsl_matrix_get(A1, j, i)));

  alphasigma = gsl_vector_alloc(6);
  for(double s=-4.0; s<=4.0; s+=0.1)
    sRepresentation(s);
  gsl_vector_free(alphasigma);

  alphasigma = gsl_vector_alloc(9);
  for(double s=-4.0; s<=4.0; s+=0.1)
    sRepresentationAsymm(s);
  gsl_vector_free(alphasigma);

  gsl_matrix_free(A1);
  gsl_matrix_free(A);
  gsl_vector_free(alpha);
}
