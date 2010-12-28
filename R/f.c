#include <gsl/gsl_vector.h>

double f(void* data, const gsl_vector* x) {
  double ans = 0.0;
  for(size_t i = 0; i < x->size; i++) {
    double xi = gsl_vector_get(x, i);
    ans -= (xi * xi)/2.0;
  }
  return(ans);
}
