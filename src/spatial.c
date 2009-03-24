#include <math.h>
#include <gsl/gsl_math.h>
#include "spatial.h"

mcmclib_spatial_lpdf* mcmclib_spatial_lpdf_alloc(gsl_vector* mu,
						 gsl_vector* rho,
						 gsl_vector* sigma,
						 gsl_vector* tausq,
						 gsl_matrix* D) {
  mcmclib_spatial_lpdf* a = (mcmclib_spatial_lpdf*) malloc(sizeof(mcmclib_spatial_lpdf));
  a->mu = mu;
  a->rho = rho;
  a->sigma = sigma;
  a->tausq = tausq;
  a->D = D;
  a->Sigma = gsl_matrix_alloc(mu->size, mu->size);
  a->norm = mcmclib_mvnorm_lpdf_alloc(mu, a->Sigma->data);
  return a;
}

void mcmclib_spatial_lpdf_free(mcmclib_spatial_lpdf* p) {
  gsl_matrix_free(p->Sigma);
  mcmclib_mvnorm_lpdf_free(p->norm);
  free(p);
}

double mcmclib_spatial_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_spatial_lpdf* p = (mcmclib_spatial_lpdf*) in_p;
  int d = p->D->size1;
  for(int i=0; i<d; i++)
    for(int j=0; j<d; j++) {
      gsl_matrix_set(p->Sigma, i, j, 0.0);
    }
  gsl_matrix_set_identity(p->Sigma);
  return mcmclib_mvnorm_lpdf_compute(p->norm, x);
}
