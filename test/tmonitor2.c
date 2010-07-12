/**Test monitor object: ECDF module*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <monitor.h>

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

int main() {
  gsl_matrix* X0 = gsl_matrix_alloc(10, 2);
  for(size_t i=0; i<X0->size1; i++) {
    gsl_vector_view X0i_v = gsl_matrix_row(X0, i);
    gsl_vector_set_all(&(X0i_v.vector), (double)i / (double)X0->size1);
  }
  mcmclib_monitor_ecdf* p = mcmclib_monitor_ecdf_alloc(X0);
  gsl_matrix_free(X0);
  /*check state right after init*/
  assert(gsl_vector_isnull(p->Fn));

  /*check state right after some updates*/
  gsl_vector* x = gsl_vector_alloc(2);
  gsl_vector_set_all(x, -0.1);
  mcmclib_monitor_ecdf_update(p, x);
  assert(gsl_vector_isnull(p->Fn));

  gsl_vector_set_all(x, 1.1);
  mcmclib_monitor_ecdf_update(p, x);
  for(size_t i=0; i< p->Fn->size; i++)
    assert(gsl_vector_get(p->Fn, i) == 0.5);

  gsl_vector_set_all(x, 0.5);
  mcmclib_monitor_ecdf_update(p, x);
  for(size_t i=0; i<6; i++)
    assert(check_dequal(gsl_vector_get(p->Fn, i), 2.0 / 3.0));
  for(size_t i=6; i<10; i++)
    assert(check_dequal(gsl_vector_get(p->Fn, i), 1.0 / 3.0));

  mcmclib_monitor_ecdf_free(p);
  gsl_vector_free(x);
  return 0;
}
