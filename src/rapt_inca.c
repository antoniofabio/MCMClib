#include "rapt_inca.h"
#include "vector_stats.h"

mcmclib_inca_rapt* mcmclib_inca_rapt_alloc(gsl_rng* r,
				      distrfun_p logdistr, void* logdistr_data,
				      gsl_vector** x, int M, int t0,
				      const gsl_matrix* sigma_whole, int K,
				      gsl_matrix** sigma_local,
				      region_fun_t which_region,
				      void* which_region_data){
  mcmclib_inca_rapt* ans = (mcmclib_inca_rapt*) malloc(sizeof(mcmclib_inca_rapt));
  int dim = x[0]->size;
  ans->M = M;
  ans->sm = (mcmclib_rapt**) malloc(M * sizeof(mcmclib_rapt*));
  for(int m=0; m<M; m++) {
    ans->sm[m] = mcmclib_rapt_alloc(r, logdistr, logdistr_data, x[m],
				    t0, sigma_whole, K, sigma_local,
				    which_region, which_region_data);
  }

  ans->t = 0;
  ans->n = gsl_vector_alloc(K);
  gsl_vector_set_all(ans->n, 0.0);
  ans->global_mean = gsl_vector_alloc(dim);
  ans->global_variance = gsl_matrix_alloc(dim, dim);
  ans->means = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  ans->variances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    ans->means[k] = gsl_vector_alloc(dim);
    ans->variances[k] = gsl_matrix_alloc(dim, dim);
  }
  return ans;
}

void mcmclib_inca_rapt_free(mcmclib_inca_rapt* p) {
  int K = p->sm[0]->K;
  for(int m=0; m<p->M; m++)
    mcmclib_rapt_free(p->sm[m]);
  free(p->sm);
  gsl_vector_free(p->n);
  for(int k=0; k<K; k++) {
    gsl_vector_free(p->means[k]);
    gsl_matrix_free(p->variances[k]);
  }
  free(p->means);
  free(p->variances);
  gsl_vector_free(p->global_mean);
  gsl_matrix_free(p->global_variance);
  free(p);
}

int mcmclib_inca_rapt_update(mcmclib_inca_rapt* p) {
  gsl_vector** means = p->means;
  gsl_vector* global_mean = p->global_mean;
  gsl_matrix** variances = p->variances;
  gsl_matrix* global_variance = p->global_variance;

  for(int m=0; m < p->M; m++) {
    mcmclib_rapt* s = p->sm[m];
    mcmclib_rapt_update(s);

    /* update common means and variances */
    mcmclib_covariance_update(global_variance, global_mean, &(p->t), s->current_x);
    int k = s->which_region_x;
    gsl_vector_set(p->n, k, gsl_vector_get(p->n, k) + 1);
    int fake = round(gsl_vector_get(p->n, k)) - 1;
    mcmclib_covariance_update(variances[k], means[k], &fake, s->current_x);

    if(s->t > s->t0)
      mcmclib_rapt_update_proposals_custom(s, variances, global_variance);
  }
  return 1;
}
