#include <gsl/gsl_blas.h>
#include <assert.h>
#include "gauss_am.h"
#include "gauss_mrw.h"

mcmclib_amh* mcmclib_gauss_am_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* logdistr_data,
				    gsl_vector* start_x,
				    const gsl_matrix* sigma_zero, int t0) {
  mcmclib_gauss_am_suff* suff = (mcmclib_gauss_am_suff*)
    malloc(sizeof(mcmclib_gauss_am_suff*));
  int d = start_x->size;
  suff->Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_memcpy(suff->Sigma_zero, sigma_zero);
  suff->t0 = t0;
  suff->sum_x = gsl_vector_alloc(d);
  gsl_vector_set_zero(suff->sum_x);
  suff->sum_xx = gsl_matrix_alloc(d, d);
  gsl_matrix_set_zero(suff->sum_xx);
  suff->sf = (2.38 * 2.38) / (double) d;

  suff->Sigma_eps = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(suff->Sigma_eps);
  gsl_matrix_scale(suff->Sigma_eps, 0.001);

  return mcmclib_amh_alloc(mcmclib_gauss_mrw_alloc(r, logdistr, logdistr_data,
						   start_x, sigma_zero),
			   suff, mcmclib_gauss_am_update_gamma);
}

void mcmclib_gauss_am_free(mcmclib_amh* p) {
  mcmclib_gauss_am_suff* s = p->suff;
  gsl_matrix_free(s->Sigma_eps);
  gsl_matrix_free(s->Sigma_zero);
  gsl_vector_free(s->sum_x);
  gsl_matrix_free(s->sum_xx);
  free(s);
  mcmclib_gauss_mrw_free(p->mh);
  mcmclib_amh_free(p);
}

void mcmclib_gauss_am_update_suff(mcmclib_gauss_am_suff* s, gsl_vector* x) {
  gsl_vector_add(s->sum_x, x);
  gsl_matrix_view x_cv = gsl_matrix_view_array(x->data, x->size, 1);
  gsl_matrix* x_cm = &(x_cv.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, x_cm, x_cm, 1.0, s->sum_xx);
}

void mcmclib_gauss_am_update_gamma(void* in_p, gsl_vector* x) {
  mcmclib_amh* p = (mcmclib_amh*) in_p;
  mcmclib_gauss_am_suff* s = p->suff;
  int t0 = s->t0;
  int t = p->n;

  mcmclib_gauss_am_update_suff(s, x);

  if(t >= t0) {
    mcmclib_gauss_mrw_gamma* g = (mcmclib_gauss_mrw_gamma*) p->mh->q->gamma;
    int d = x->size;
    gsl_vector* mean = gsl_vector_alloc(d);
    gsl_vector_memcpy(mean, s->sum_x);
    gsl_vector_scale(mean, 1.0 / (double) t);
    gsl_matrix_view mean_cv = gsl_matrix_view_array(mean->data, d, 1);
    gsl_matrix* mean_cm = &(mean_cv.matrix);

    gsl_matrix_memcpy(g->Sigma, s->sum_xx);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, mean_cm, mean_cm, 1.0, g->Sigma);
    gsl_matrix_scale(g->Sigma, s->sf);
    gsl_vector_free(mean);
  }
}

void mcmclib_gauss_am_set_sf(mcmclib_amh* p, double sf) {
  assert(sf > 0.0);
  ((mcmclib_gauss_am_suff*) (p->suff))->sf = sf;
}

void mcmclib_gauss_am_set_eps(mcmclib_amh* p, double eps) {
  assert(eps > 0.0);
  mcmclib_gauss_am_suff* s = (mcmclib_gauss_am_suff*) p->suff;
  gsl_matrix_set_identity(s->Sigma_eps);
  gsl_matrix_scale(s->Sigma_eps, eps);
}
