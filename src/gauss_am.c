#include "gauss_am.h"
#include "vector_stats.h"
#include "gauss_mrw.h"

mcmclib_gauss_am* mcmclib_gauss_am_alloc(gsl_rng* r,
					 distrfun_p logdistr, void* logdistr_data,
					 gsl_vector* start_x,
					 const gsl_matrix* sigma_zero, int t0) {
  mcmclib_gauss_am* a = (mcmclib_gauss_am*) malloc(sizeof(mcmclib_gauss_am));
  a->mrw = mcmclib_gauss_mrw_alloc(r, logdistr, logdistr_data, start_x, sigma_zero);
  a->amh = mcmclib_amh_alloc(a->mrw->mh, mcmclib_gauss_am_update_gamma, a);
  int d = sigma_zero->size1;

  a->sigma_zero = gsl_matrix_alloc(d, d);
  a->t0 = t0;
  a->mean = gsl_vector_alloc(d);
  a->cov = gsl_matrix_alloc(d, d);
  gsl_matrix_memcpy(a->sigma_zero, sigma_zero);
  a->sf = (2.38 * 2.38) / d;

  gsl_vector_set_zero(a->mean);
  gsl_matrix_set_zero(a->cov);
  return a;
}

void mcmclib_gauss_am_free(mcmclib_gauss_am* p) {
  mcmclib_gauss_mrw_free(p->mrw);
  gsl_matrix_free(p->sigma_zero);
  gsl_vector_free(p->mean);
  gsl_matrix_free(p->cov);
  free(p);
}

int mcmclib_gauss_am_update(mcmclib_gauss_am* p) {
  return mcmclib_amh_update(p->amh);
}

void mcmclib_gauss_am_update_gamma(void* in_p) {
  mcmclib_gauss_am* p = (mcmclib_gauss_am*) in_p;
  int t0 = p->t0;
  int t = p->amh->n;
  mcmclib_gauss_mrw* mrw = p->mrw;

  mcmclib_covariance_update(p->cov, p->mean, &t, p->amh->mh->x);
  if(t >= t0) {
    gsl_matrix_memcpy(mrw->sigma_prop, p->cov);
    gsl_matrix_scale(mrw->sigma_prop, p->sf);
  }
}
