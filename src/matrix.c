#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include "matrix.h"

int mcmclib_cholesky_decomp(gsl_matrix* A) {
  gsl_error_handler_t *hnd = gsl_set_error_handler_off();
  int status = gsl_linalg_cholesky_decomp(A);
  gsl_set_error_handler(hnd);
  if(status != GSL_SUCCESS)
    return status;
  return GSL_SUCCESS;
}

int mcmclib_cholesky_inverse(gsl_matrix* A) {
  int status = mcmclib_cholesky_decomp(A);
  if(status)
    return status;
  gsl_linalg_cholesky_invert(A);
  return GSL_SUCCESS;
}

double mcmclib_matrix_logtrace(const gsl_matrix* A) {
  double ans = 0.0;
  for(int i=0; i < A->size1; i++)
    ans += log(gsl_matrix_get(A, i, i));
  return ans;
}

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

int mcmclib_vector_is_sorted_desc(gsl_vector* v) {
  double m = gsl_vector_get(v, 0);
  for(int i=1; i<v->size; i++) {
    double n = gsl_vector_get(v, i);
    if(n > m)
      return 0;
    m = n;
  }
  return 1;
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
