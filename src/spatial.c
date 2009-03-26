#include <assert.h>
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
  assert(mu->size == D->size1);
  a->rho = rho;
  assert(rho->size == 1);
  a->sigma = sigma;
  assert(sigma->size == 1);
  a->tausq = tausq;
  assert(tausq->size == 1);
  a->D = D;
  assert(D->size1 == D->size2);
  a->Sigma = gsl_matrix_alloc(mu->size, mu->size);
  a->norm = mcmclib_mvnorm_lpdf_alloc(mu, a->Sigma->data);
  return a;
}

void mcmclib_spatial_lpdf_free(mcmclib_spatial_lpdf* p) {
  gsl_matrix_free(p->Sigma);
  mcmclib_mvnorm_lpdf_free(p->norm);
  free(p);
}

double mcmclib_spatial_cov_exponential(double d, double rho, double sigma, double tausq) {
  double a = (sigma - tausq) * (1.0 - exp(- d / rho));
  if (d>0)
    a += tausq;
  return a;
}

void mcmclib_spatial_distances(gsl_matrix* D, gsl_matrix* xy) {
  int n = xy->size1;
  int d = xy->size2;
  gsl_matrix_set_all(D, 0.0);
  gsl_vector* x = gsl_vector_alloc(d);
  for(int i=0; i<n; i++) for(int j=i+1; j<n; j++) {
      gsl_vector_view ri_v = gsl_matrix_row(xy, i);
      gsl_vector* ri = &(ri_v.vector);
      gsl_vector_memcpy(x, ri);
      gsl_vector_view rj_v = gsl_matrix_row(xy, j);
      gsl_vector* rj = &(rj_v.vector);
      gsl_vector_sub(x, rj);
      gsl_vector_mul(x, x);
      double dist = 0.0;
      for(int k = 0; k<d; k++)
	dist += gsl_vector_get(x, k);
      dist = sqrt(dist);
      gsl_matrix_set(D, i, j, dist);
      gsl_matrix_set(D, j, i, dist);
  }
  gsl_vector_free(x);
}

void mcmclib_spatial_set_xy(mcmclib_spatial_lpdf* p, gsl_matrix* xy) {
  mcmclib_spatial_distances(p->D, xy);
}

double mcmclib_spatial_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_spatial_lpdf* p = (mcmclib_spatial_lpdf*) in_p;
  int n = x->size;
  gsl_matrix* D = p->D;
  assert(n == D->size1);
  double rho = gsl_vector_get(p->rho, 0);
  double sigma = gsl_vector_get(p->sigma, 0);
  double tausq = gsl_vector_get(p->tausq, 0);
  for(int i=0; i<n; i++)
    gsl_matrix_set(p->Sigma, i, i, sigma);
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++) {
      double gammaij = mcmclib_spatial_cov_exponential(gsl_matrix_get(D, i, j), rho, sigma, tausq);
      gsl_matrix_set(p->Sigma, i, j, sigma - gammaij);
      gsl_matrix_set(p->Sigma, j, i, sigma - gammaij);
    }
  return mcmclib_mvnorm_lpdf_compute(p->norm, x);
}
