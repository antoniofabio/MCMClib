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

void mcmclib_matrix_printf(gsl_matrix* A) {
  int n = A->size1;
  int p = A->size2;
  for(int i=0; i<n; i++) {
    for(int j=0; j<p; j++)
      printf("%.3f, ", gsl_matrix_get(A, i, j));
    printf("\n");
  }
}
