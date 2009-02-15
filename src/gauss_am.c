#include "gauss_am.h"
#include "vector_stats.h"
#include "gauss_mrw.h"

mcmclib_gauss_am* mcmclib_gauss_am_alloc(gsl_rng* r,
					 distrfun_p logdistr, void* logdistr_data,
					 gsl_vector* start_x,
					 const gsl_matrix* sigma_zero, int t0) {
  mcmclib_gauss_am* ans = (mcmclib_gauss_am*) malloc(sizeof(mcmclib_gauss_am));
  mcmclib_gauss_mrw* mrw = mcmclib_gauss_mrw_alloc(r, logdistr, logdistr_data, start_x, sigma_zero);
  ans->mrw = mrw;
  int d = sigma_zero->size1;
  ans->r = r;
  ans->logdistr = logdistr;
  ans->logdistr_data = logdistr_data;
  ans->current_x = mrw->mh->x;
  ans->old = mrw->mh->x_old;

  ans->sigma_zero = gsl_matrix_alloc(d, d);
  ans->t0 = t0;
  ans->mean = gsl_vector_alloc(d);
  ans->cov = gsl_matrix_alloc(d, d);
  ans->t = 0;
  gsl_matrix_memcpy(ans->sigma_zero, sigma_zero);
  ans->sf = (2.38 * 2.38) / d;

  gsl_vector_set_zero(ans->mean);
  gsl_matrix_set_zero(ans->cov);
  return ans;
}

void mcmclib_gauss_am_free(mcmclib_gauss_am* p) {
	mcmclib_gauss_mrw_free(p->mrw);
	gsl_matrix_free(p->sigma_zero);
	gsl_vector_free(p->mean);
	gsl_matrix_free(p->cov);
	free(p);
}

int mcmclib_gauss_am_update(mcmclib_gauss_am* p) {
	int t0 = p->t0;
	int *t = &(p->t);
	mcmclib_gauss_mrw* mrw = p->mrw;
	gsl_matrix* cov = p->cov;

	/*sample a new value for 'x', with mean 'old' and assigned covariance matrix*/
	mcmclib_gauss_mrw_update(mrw);

	/*adapt extra parameters*/
	mcmclib_covariance_update(p->cov, p->mean, t, p->current_x);
	if(((*t) + 1) >= t0) {
		gsl_matrix_memcpy(mrw->sigma_prop, cov);
		gsl_matrix_scale(mrw->sigma_prop, p->sf);
	}

	return 0;
}
