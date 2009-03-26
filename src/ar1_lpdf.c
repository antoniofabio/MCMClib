#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "ar1_lpdf.h"

mcmclib_ar1_lpdf* mcmclib_spatial_lpdf_alloc(gsl_vector* mu,
					     gsl_vector* sigma,
					     gsl_vector* phi) {
  mcmclib_ar1_lpdf* a = (mcmclib_ar1_lpdf*) malloc(sizeof(mcmclib_ar1_lpdf));
  a->mu = mu;
  a->sigma = sigma;
  assert(sigma->size==1);
  a->phi = phi;
  assert(phi->size==1);
  a->norm = NULL; //TODO
  return a;
}

void mcmclib_ar1_lpdf_free(mcmclib_ar1_lpdf* p) {
  free(p);
}

double mcmclib_ar1_lpdf_compute(void* in_p, gsl_vector* x) {
  return 0.0;
}
