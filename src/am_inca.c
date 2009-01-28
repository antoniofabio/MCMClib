#include "am_inca.h"
#include "vector_stats.h"

mcmclib_am_inca* mcmclib_am_inca_alloc(gsl_rng* r,
				       distrfun_p logdistr, void* data,
				       gsl_vector** x, int M, int t0,
				       const gsl_matrix* sigma_prop) {
  mcmclib_am_inca* ans = (mcmclib_am_inca*) malloc(sizeof(mcmclib_am_inca));
  int dim = x[0]->size;
  ans->M = M;
  ans->sm = (mcmclib_gauss_mrw**) malloc(M * sizeof(mcmclib_gauss_mrw*));
  for(int m=0; m<M; m++) {
    ans->sm[m] = mcmclib_gauss_mrw_alloc(r, logdistr, data, x[m], sigma_prop);
  }

  ans->t0 = t0;
  ans->t = 0;
  ans->global_mean = gsl_vector_alloc(dim);
  ans->global_variance = gsl_matrix_alloc(dim, dim);
  ans->Sigma_eps = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_identity(ans->Sigma_eps);
  gsl_matrix_scale(ans->Sigma_eps, 0.001);
  return ans;
}

void mcmclib_am_inca_free(mcmclib_am_inca* p) {
  for(int m=0; m<p->M; m++)
    mcmclib_gauss_mrw_free(p->sm[m]);
  free(p->sm);
  gsl_vector_free(p->global_mean);
  gsl_matrix_free(p->global_variance);
  free(p);
}

int mcmclib_am_inca_update(mcmclib_am_inca* p) {
  double dim = (double) p->sm[0]->current_x->size;
  gsl_vector* global_mean = p->global_mean;
  gsl_matrix* global_variance = p->global_variance;
  int M = p->M;

  for(int m=0; m < M; m++) {
    mcmclib_gauss_mrw* s = p->sm[m];
    mcmclib_gauss_mrw_update(s);

    /* update common means and variances */
    mcmclib_covariance_update(global_variance, global_mean, &(p->t), s->current_x);

    if((p->t/M) > p->t0) {
      gsl_matrix_memcpy(s->sigma_prop, global_variance);
      gsl_matrix_add(s->sigma_prop, p->Sigma_eps);
      gsl_matrix_scale(s->sigma_prop, 2.38 * 2.38 / dim);
    }
  }
  return 1;
}
