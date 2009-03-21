#include <gsl/gsl_math.h>
#include "rapt.h"
#include "mvnorm.h"
#include "vector_stats.h"

mcmclib_rapt_suff* mcmclib_rapt_suff_alloc(int t0, int K, int dim) {
  mcmclib_rapt_suff* a = (mcmclib_rapt_suff*) malloc(sizeof(mcmclib_rapt_suff));
  a->t0 = t0;

  a->means = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  a->variances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    a->means[k] = gsl_vector_alloc(dim);
    a->variances[k] = gsl_matrix_alloc(dim, dim);
  }
  a->global_mean = gsl_vector_alloc(dim);
  gsl_vector_set_all(a->global_mean, 0.0);
  a->global_variance = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_all(a->global_variance, 0.0);
  a->n = gsl_vector_alloc(K);
  gsl_vector_set_all(a->n, 0.0);
  a->Sigma_eps = gsl_matrix_alloc(dim, dim);

  mcmclib_rapt_suff_set_correction_factor(a, 0.001);
  double sf = 2.38 * 2.38 / ((double) dim);
  mcmclib_rapt_suff_set_scaling_factors(a, sf, sf);

  return a;
}

void mcmclib_rapt_suff_free(mcmclib_rapt_suff* p) {
  int K = p->n->size;
  gsl_vector_free(p->n);
  for(int k=0; k< K; k++) {
    gsl_vector_free(p->means[k]);
    gsl_matrix_free(p->variances[k]);
  }
  gsl_vector_free(p->global_mean);
  gsl_matrix_free(p->global_variance);
  free(p->means);
  free(p->variances);

  free(p);
}

void mcmclib_rapt_update(void* p);

mcmclib_amh* mcmclib_rapt_alloc(gsl_rng* r,
				distrfun_p logdistr, void* logdistr_data,
				gsl_vector* x,
				int t0,
				const gsl_matrix* sigma_whole,
				int K,
				gsl_matrix** sigma_local,
				region_fun_t which_region,
				void* which_region_data) {
  int dim = x->size;

  mcmclib_mh_q* q = mcmclib_rapt_q_alloc(r, logdistr, logdistr_data,
					 sigma_whole, K, sigma_local,
					 which_region, which_region_data);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, logdistr, logdistr_data, q, x);
  mcmclib_rapt_suff* suff = mcmclib_rapt_suff_alloc(t0, K, dim);

  return mcmclib_amh_alloc(mh, suff, mcmclib_rapt_update);
}

void mcmclib_rapt_free(mcmclib_amh* p) {
  mcmclib_rapt_q_free(p->mh->q);
  mcmclib_mh_free(p->mh);
  mcmclib_rapt_suff_free((mcmclib_rapt_suff*) p->suff);
  mcmclib_amh_free(p);
}

void mcmclib_rapt_update_proposals_custom(mcmclib_amh* p,
					  gsl_matrix** variances,
					  gsl_matrix* global_variance) {
  mcmclib_rapt_gamma* g = (mcmclib_rapt_gamma*) p->mh->q->gamma;
  mcmclib_rapt_suff* s = (mcmclib_rapt_suff*) p->suff;
  mcmclib_rapt_q_update_proposals_custom(g, variances, global_variance,
					 s->Sigma_eps,
					 s->scaling_factor_local,
					 s->scaling_factor_global);
}

void mcmclib_rapt_update_proposals(mcmclib_amh* p) {
  mcmclib_rapt_suff* s = (mcmclib_rapt_suff*) p->suff;
  if((p->n) <= s->t0)
    return;
  mcmclib_rapt_update_proposals_custom(p, s->variances, s->global_variance);
}

void mcmclib_rapt_update_suff(mcmclib_amh* p) {
  mcmclib_rapt_suff* s = (mcmclib_rapt_suff*) p->suff;
  mcmclib_mh* mh = p->mh;
  mcmclib_rapt_gamma* g = (mcmclib_rapt_gamma*) mh->q->gamma;
  gsl_vector* x = mh->x;
  int k = g->which_region_x = g->which_region(x, g->which_region_data);
  int fake_n = gsl_vector_get(s->n, k);
  mcmclib_covariance_update(s->variances[k], s->means[k], &fake_n, x);
  fake_n = p->n - 1;
  mcmclib_covariance_update(s->global_variance, s->global_mean, &fake_n, x);
  gsl_vector_set(s->n, g->which_region_x, gsl_vector_get(s->n, g->which_region_x) + 1);
}

void mcmclib_rapt_update(void* in_p) {
  mcmclib_amh* p = (mcmclib_amh*) in_p;
  mcmclib_rapt_update_suff(p);
  mcmclib_rapt_update_proposals(p);
}

void mcmclib_rapt_suff_set_correction_factor(mcmclib_rapt_suff* p, double eps) {
  gsl_matrix_set_identity(p->Sigma_eps);
  gsl_matrix_scale(p->Sigma_eps, eps);
}

void mcmclib_rapt_suff_set_scaling_factors(mcmclib_rapt_suff* p, double local, double global) {
  p->scaling_factor_local = local;
  p->scaling_factor_global = global;
}
