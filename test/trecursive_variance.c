/**Test recursive variance function*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <vector_stats.h>

#define DIM 2
#define N 3

int main(int argc, char** argv) {
  gsl_vector* x[N];
  for(int n=0; n<N; n++) {
    x[n] = gsl_vector_alloc(DIM);
    for(int d=0; d<DIM; d++)
      gsl_vector_set(x[n], d, n*N - N + n*d);
  }

  gsl_vector* m = gsl_vector_alloc(DIM);
  gsl_matrix* V = gsl_matrix_alloc(DIM, DIM);

  int n1=0;
  for(int n=0; n<N; n++) {
    mcmclib_covariance_update(V, m, &n1, x[n]);
  }

  gsl_matrix* Vt = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix* X = gsl_matrix_alloc(N, DIM);
  for(int n=0; n<N; n++) {
    gsl_vector_view Xn = gsl_matrix_row(X, n);
    gsl_vector_memcpy(&(Xn.vector), x[n]);
  }
  mcmclib_matrix_covariance(X, Vt);

  /*check results*/
  assert(n1 == N);
  for(int d1=0; d1<DIM; d1++)
    for(int d2=0; d2<DIM; d2++)
      assert(gsl_matrix_get(V, d1, d2) == gsl_matrix_get(Vt, d1, d2));

  /*free memory*/
  for(int n=0; n<N; n++)
    gsl_vector_free(x[n]);
  gsl_matrix_free(V);
  gsl_matrix_free(Vt);
  gsl_matrix_free(X);
  gsl_vector_free(m);

  return 0;
}
