/**Test recursive variance function*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <vector_stats.h>

#define DIM 2
#define N 3

int main() {
  gsl_vector* x[N];
  for(size_t n=0; n<N; n++) {
    x[n] = gsl_vector_alloc(DIM);
    for(size_t d=0; d<DIM; d++)
      gsl_vector_set(x[n], d, (double) (n*N) - (double) N + (double) (n*d));
  }

  gsl_vector* m = gsl_vector_alloc(DIM);
  gsl_matrix* V = gsl_matrix_alloc(DIM, DIM);

  size_t n1=0;
  for(size_t n=0; n<N; n++) {
    mcmclib_covariance_update(V, m, &n1, x[n]);
  }

  gsl_matrix* Vt = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  for(size_t n=0; n<N; n++) {
    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x[n]);
  }
  mcmclib_matrix_covariance(X, Vt);

  /*check results*/
  assert(n1 == N);
  for(size_t d1=0; d1<DIM; d1++)
    for(size_t d2=0; d2<DIM; d2++)
      assert(gsl_matrix_get(V, d1, d2) == gsl_matrix_get(Vt, d1, d2));

  /*free memory*/
  for(size_t n=0; n<N; n++)
    gsl_vector_free(x[n]);
  gsl_matrix_free(V);
  gsl_matrix_free(Vt);
  gsl_matrix_free(X);
  gsl_vector_free(m);

  return 0;
}
