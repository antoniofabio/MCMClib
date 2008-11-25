/** \file
\brief Compute distance from a target distribution

Target distrib. given as a sample of obs, puts output to stdout
 */

#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/**Compute Fn in point 'z', basing on sample 'X'*/
double Fn(gsl_vector* z, gsl_matrix* X, gsl_vector* workspace) {
  int N = X->size1;
  double ans = 0.0;
  for(int n=0; n < N; n++) {
    gsl_vector_memcpy(workspace, z);
    gsl_vector_view xv = gsl_matrix_row(X, n);
    gsl_vector* x = &(xv.vector);
    gsl_vector_sub(workspace, x);
    ans += !gsl_vector_ispos(workspace);
  }
  return ans / (double) N;
}

/**Compute Fn for each row in 'X0'*/
gsl_vector* Fn_vector(gsl_matrix* X, gsl_matrix* X0) {
  assert(X->size2 == X0->size2);
  gsl_vector* a = gsl_vector_alloc(X0->size1);
  gsl_vector* workspace = gsl_vector_alloc(X0->size2);
  for(int n=0; n < X0->size1; n++) {
    gsl_vector_view zv = gsl_matrix_row(X0, n);
    gsl_vector_set(a, n, Fn(&(zv.vector), X, workspace));
  }
  gsl_vector_free(workspace);
  return a;
}

gsl_vector* vmap(gsl_vector* x, double (*op)(double)) {
  for(int i=0; i < x->size; i++)
    gsl_vector_set(x, i, op(gsl_vector_get(x, i)));
  return x;
}

/**Squares all elements of 'x', in place*/
double square(double x) {return x*x;}
gsl_vector* vector_square(gsl_vector* x) {
  return vmap(x, square);
}

int main(int argc, char** argv) {
  if(argc < 8) {
    printf("usage: %s dim X_filename X_N Y_filename Y_N X0_filename X0_N\n", argv[0]);
    return 1;
  }
  gsl_set_error_handler_off();
  int dim, X_N;
  sscanf(argv[1], "%d", &dim);
  FILE* f = fopen(argv[2], "r");
  sscanf(argv[3], "%d", &X_N);
  gsl_matrix* X = gsl_matrix_alloc(X_N, dim);
  gsl_matrix_fscanf(f, X);
  fclose(f);
  int Y_N;
  f = fopen(argv[4], "r");
  sscanf(argv[5], "%d", &Y_N);
  gsl_matrix* Y = gsl_matrix_alloc(Y_N, dim);
  gsl_matrix_fscanf(f, Y);
  fclose(f);
  int X0_N;
  f = fopen(argv[6], "r");
  sscanf(argv[7], "%d", &X0_N);
  gsl_matrix* X0 = gsl_matrix_alloc(X0_N, dim);
  gsl_matrix_fscanf(f, X0);
  fclose(f);

  /*main loop*/
  gsl_vector* Fn_Y = Fn_vector(Y, X0);
  for(int n=0; n < X_N; n++) {
    gsl_matrix_view XXn = gsl_matrix_submatrix(X, 0, 0, n+1, dim);
    gsl_vector* Fn_XXn = Fn_vector(&(XXn.matrix), X0);
    gsl_vector_sub(Fn_XXn, Fn_Y);
    vector_square(Fn_XXn);
    printf("%f\n", gsl_stats_mean(Fn_XXn->data, Fn_XXn->stride, X0_N));
    gsl_vector_free(Fn_XXn);
  }
  gsl_vector_free(Fn_Y);
  /*end main loop*/

  gsl_matrix_free(X0);
  gsl_matrix_free(Y);
  gsl_matrix_free(X);
}
