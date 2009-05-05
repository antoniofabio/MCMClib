#include <gsl/gsl_linalg.h>
#include "matrix.h"

void mcmclib_matrix_inverse(gsl_matrix* A) {
  gsl_permutation* p = gsl_permutation_alloc(A->size1);
  gsl_matrix* A1 = gsl_matrix_alloc(A->size1, A->size1);
  int tmp=0;
  gsl_matrix_memcpy(A1, A);
  gsl_linalg_LU_decomp(A1, p, &tmp);
  gsl_linalg_LU_invert(A1, p, A);
  gsl_matrix_free(A1);
  gsl_permutation_free(p);
}

void mcmclib_vector_printf(gsl_vector* v) {
  int n = v->size;
  printf("%.3f", gsl_vector_get(v, 0));
  for(int i=1; i<n; i++) {
    printf(", %.3f", gsl_vector_get(v, i));
  }
  printf("\n");
}

void mcmclib_matrix_printf(gsl_matrix* A) {
  int n = A->size1;
  for(int i=0; i<n; i++) {
    gsl_vector_view row = gsl_matrix_row(A, i);
    mcmclib_vector_printf(&row.vector);
  }
}
