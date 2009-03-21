#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include "mixem_online.h"
#include "vector_stats.h"
#include "mvnorm.h"

mcmclib_mixem_online* mcmclib_mixem_online_alloc(gsl_vector** mu,
						 gsl_matrix** Sigma,
						 gsl_vector* beta,
						 double eta_eps,
						 int n0) {
  mcmclib_mixem_online* a = (mcmclib_mixem_online*) malloc(sizeof(mcmclib_mixem_online));
  a->beta = beta;
  a->mu = mu;
  a->Sigma = Sigma;
  a->gamma = mcmclib_mixolem_suff_makeref(beta, mu, Sigma);
  assert((eta_eps > 0) && (eta_eps <= 0.5));
  a->eta_eps = eta_eps;
  a->n0 = n0;
  a->n = 0;

  int K = beta->size;
  int d = mu[0]->size;
  a->s = mcmclib_mixolem_suff_alloc(K, d);
  a->si = mcmclib_mixolem_suff_alloc(K, d);

  a->pi_k = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++) {
    a->pi_k[k] = mcmclib_mvnorm_lpdf_alloc(mu[k], Sigma[k]->data);
  }

  a->mu_global = gsl_vector_alloc(d);
  gsl_vector_set_all(a->mu_global, 0.0);
  a->Sigma_global = gsl_matrix_alloc(d, d);
  gsl_matrix_set_all(a->Sigma_global, 0.0);
  return a;
}

void mcmclib_mixem_online_free(mcmclib_mixem_online* p) {
  for(int k=0; k < p->gamma->delta->size; k++) {
    mcmclib_mvnorm_lpdf_free(p->pi_k[k]);
  }
  free(p->gamma);
  mcmclib_mixolem_suff_free(p->s);
  mcmclib_mixolem_suff_free(p->si);
  gsl_vector_free(p->mu_global);
  gsl_matrix_free(p->Sigma_global);
  free(p->pi_k);
  free(p);
}

/*update sufficient statistic*/
static void mixem_compute_si(mcmclib_mixem_online* p, gsl_vector* y) {
  int K = p->si->delta->size;
  double delta_sum = 0.0;
  for(int k=0; k < K; k++) {
    double delta_k = 0.0;
    delta_k = exp(mcmclib_mvnorm_lpdf_compute(p->pi_k[k], y) +
      log(gsl_vector_get(p->gamma->delta, k)));
    gsl_vector_set(p->si->delta, k, delta_k);
    delta_sum += delta_k;
  }
  gsl_vector_scale(p->si->delta, 1.0 / delta_sum);

  for(int k=0; k< K; k++) {
    double delta_k = gsl_vector_get(p->si->delta, k);
    gsl_vector_memcpy(p->si->delta_x[k], y);
    gsl_vector_scale(p->si->delta_x[k], delta_k);

    gsl_matrix_view xmv = gsl_matrix_view_array(y->data, y->size, 1);
    gsl_matrix* xm = &(xmv.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, xm, xm,
		   0.0, p->si->delta_xx[k]);
    gsl_matrix_scale(p->si->delta_xx[k], delta_k);
  }
}
void mcmclib_mixem_online_update_s(mcmclib_mixem_online* p, gsl_vector* y) {
  mixem_compute_si(p, y);
  double eta_n = pow((double) p->n, -0.5 - p->eta_eps);
  mcmclib_mixolem_suff_scale(p->si, eta_n);
  mcmclib_mixolem_suff_scale(p->s, 1.0 - eta_n);
  mcmclib_mixolem_suff_add(p->s, p->si);
}

void mcmclib_mixem_online_update_gamma(mcmclib_mixolem_suff* gamma,
				       mcmclib_mixolem_suff* s) {
  mcmclib_mixolem_suff_memcpy(gamma, s);
  int K = gamma->delta->size;
  for(int k=0; k < K; k++) {
    double betak = 1.0 / gsl_vector_get(gamma->delta, k);
    gsl_vector_scale(gamma->delta_x[k], betak);

    gsl_matrix* Sigmak = gamma->delta_xx[k];
    gsl_matrix_scale(Sigmak, betak);
    gsl_vector* muk = gamma->delta_x[k];
    gsl_matrix_view mucv = gsl_matrix_view_array(muk->data, muk->size, 1);
    gsl_matrix* muc = &(mucv.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, muc, muc, 1.0, Sigmak);
  }
}

void mcmclib_mixem_online_update(mcmclib_mixem_online* p, gsl_vector* y) {
  int n = p->n;
  mcmclib_covariance_update(p->Sigma_global, p->mu_global, &n, y);
  (p->n)++;
  mcmclib_mixem_online_update_s(p, y);
  if(p->n <= p->n0)
    return;
  mcmclib_mixem_online_update_gamma(p->gamma, p->s);
}
