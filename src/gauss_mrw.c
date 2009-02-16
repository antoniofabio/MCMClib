#include <gsl/gsl_math.h>
#include "gauss_mrw.h"
#include "mvnorm.h"

double mcmclib_gauss_mrw_qd(void* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

void mcmclib_gauss_mrw_sample(void* in_p, gsl_vector* x) {
  mcmclib_gauss_mrw* p = (mcmclib_gauss_mrw*) in_p;
  mcmclib_mvnorm(p->mh->r, p->sigma_prop, x);
}

mcmclib_gauss_mrw* mcmclib_gauss_mrw_alloc(gsl_rng* r,
					   distrfun_p logdistr, void* logdistr_data,
					   gsl_vector* start_x, const gsl_matrix* sigma_prop) {
  mcmclib_gauss_mrw* a = (mcmclib_gauss_mrw*) malloc(sizeof(mcmclib_gauss_mrw));
  int dim = start_x->size;
  a->sigma_prop = gsl_matrix_alloc(dim, dim);
  gsl_matrix_memcpy(a->sigma_prop, sigma_prop);

  a->mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, start_x,
			   mcmclib_gauss_mrw_qd, NULL,
			   mcmclib_gauss_mrw_sample, a);
  return a;
}

void mcmclib_gauss_mrw_free(mcmclib_gauss_mrw* p) {
  mcmclib_mh_free(p->mh);
  gsl_matrix_free(p->sigma_prop);
  free(p);
}

int mcmclib_gauss_mrw_update(mcmclib_gauss_mrw* p) {
  return mcmclib_mh_update(p->mh);
}
