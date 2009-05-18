/**Test monitor object*/
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

int main(int argc, char** argv) {
  gsl_vector* x = gsl_vector_alloc(1);
  gsl_vector_set(x, 0, 2.0);
  mcmclib_monitor* p = mcmclib_monitor_alloc(x);
  gsl_vector* tmp = gsl_vector_alloc(1);
  mcmclib_monitor_update(p);

  /*check state right after init*/
  mcmclib_monitor_get_means(p, tmp);
  assert(v0(tmp) == x0);
  mcmclib_monitor_get_vars(p, tmp);
  assert(check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_ar(p, tmp);
  assert(!gsl_finite(v0(tmp)));
  mcmclib_monitor_get_msjd(p, tmp);
  assert(!gsl_finite(v0(tmp)));

  /*check state after 1 rejection*/
  mcmclib_monitor_update(p);

  mcmclib_monitor_get_means(p, tmp);
  assert(check_dequal(v0(tmp), x0));
  mcmclib_monitor_get_vars(p, tmp);
  assert(check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_ar(p, tmp);
  assert(check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_msjd(p, tmp);
  assert(check_dequal(v0(tmp), 0.0));

  /*check state after 1 rejection and 1 acceptance*/
  gsl_vector_set(x, 0, -4.0);
  mcmclib_monitor_update(p);

  mcmclib_monitor_get_means(p, tmp);
  assert(check_dequal(v0(tmp), 0.0));
  mcmclib_monitor_get_vars(p, tmp);
  assert(check_dequal(v0(tmp), 8.0));
  mcmclib_monitor_get_ar(p, tmp);
  assert(check_dequal(v0(tmp), 0.5));
  mcmclib_monitor_get_msjd(p, tmp);
  assert(check_dequal(v0(tmp), 18.0)); /*(2+4)^2 / 2*/

  mcmclib_monitor_free(p);
  gsl_vector_free(tmp);
  gsl_vector_free(x);
  return 0;
}
