#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include "mixem_online.h"
#include "mvnorm.h"

mcmclib_mixem_online* mcmclib_mixem_online_alloc(gsl_vector** mu,
						 gsl_matrix** Sigma,
						 gsl_vector* beta,
						 double eta_eps,
						 int n0) {
  mcmclib_mixem_online* a = (mcmclib_mixem_online*) malloc(sizeof(mcmclib_mixem_online));
  a->mu = mu;
  a->Sigma = Sigma;
  a->beta = beta;
  assert((eta_eps >= 0) && (eta_eps <= 0.5));
  a->eta_eps = eta_eps;
  a->n0 = n0;
  a->n = 0;

  int K = beta->size;
  int d = mu[0]->size;
  a->delta = gsl_vector_alloc(K);
  gsl_vector_set_all(a->delta, 0.0);
  a->delta_x = (gsl_vector**) malloc(sizeof(gsl_vector*) * K);
  a->delta_xx = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * K);
  for(int k=0; k<K; k++) {
    a->delta_x[k] = gsl_vector_alloc(d);
    gsl_vector_set_all(a->delta_x[k], 0.0);
    a->delta_xx[k] = gsl_matrix_alloc(d, d);
    gsl_matrix_set_all(a->delta_xx[k], 0.0);
  }
  a->deltai = gsl_vector_alloc(K);
  a->delta_xi = (gsl_vector**) malloc(sizeof(gsl_vector*) * K);
  a->delta_xxi = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * K);
  for(int k=0; k<K; k++) {
    a->delta_xi[k] = gsl_vector_alloc(d);
    a->delta_xxi[k] = gsl_matrix_alloc(d, d);
  }

  a->pi_k = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++) {
    a->pi_k[k] = mcmclib_mvnorm_lpdf_alloc(mu[k], Sigma[k]->data);
  }
  return a;
}

void mcmclib_mixem_online_free(mcmclib_mixem_online* p) {
  for(int k=0; k < p->beta->size; k++) {
    mcmclib_mvnorm_lpdf_free(p->pi_k[k]);
    gsl_vector_free(p->delta_x[k]);
    gsl_matrix_free(p->delta_xx[k]);
    gsl_vector_free(p->delta_xi[k]);
    gsl_matrix_free(p->delta_xxi[k]);
  }
  gsl_vector_free(p->delta);
  gsl_vector_free(p->deltai);
  free(p->delta_x);
  free(p->delta_xx);
  free(p->delta_xi);
  free(p->delta_xxi);
  free(p->pi_k);
  free(p);
}

/*update sufficient statistic*/
static void mixem_update_s(mcmclib_mixem_online* p, gsl_vector* y) {
  double delta_sum = 0.0;
  for(int k=0; k < p->beta->size; k++) {
    double delta_k = 0.0;
    delta_k = exp(mcmclib_mvnorm_lpdf_compute(p->pi_k[k], y) +
      log(gsl_vector_get(p->beta, k)));
    gsl_vector_set(p->deltai, k, delta_k);
    delta_sum += delta_k;
  }
  gsl_vector_scale(p->deltai, 1.0 / delta_sum);

  for(int k=0; k< p->beta->size; k++) {
    double delta_k = gsl_vector_get(p->deltai, k);
    gsl_vector_memcpy(p->delta_xi[k], y);
    gsl_vector_scale(p->delta_xi[k], delta_k);

    gsl_matrix_view xmv = gsl_matrix_view_array(y->data, y->size, 1);
    gsl_matrix* xm = &(xmv.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, xm, xm,
		   0.0, p->delta_xxi[k]);
    gsl_matrix_scale(p->delta_xxi[k], delta_k);
  }

  double eta_n = pow((double) p->n, -0.5 - p->eta_eps);

  gsl_vector_scale(p->delta, 1.0 - eta_n);
  gsl_vector_scale(p->deltai, eta_n);
  gsl_vector_add(p->delta, p->deltai);
  for(int k=0; k < p->beta->size; k++) {
    gsl_vector_scale(p->delta_x[k], 1.0 - eta_n);
    gsl_vector_scale(p->delta_xi[k], eta_n);
    gsl_vector_add(p->delta_x[k], p->delta_xi[k]);

    gsl_matrix_scale(p->delta_xx[k], 1.0 - eta_n);
    gsl_matrix_scale(p->delta_xxi[k], eta_n);
    gsl_matrix_add(p->delta_xx[k], p->delta_xxi[k]);
  }
}

/*update ML estimate*/
static void mixem_update_gamma(mcmclib_mixem_online* p) {
  if(p->n <= p->n0)
    return;
  gsl_vector_memcpy(p->beta, p->delta);
  for(int k=0; k < p->beta->size; k++) {
    double betak = 1.0 / gsl_vector_get(p->beta, k);
    gsl_vector_memcpy(p->mu[k], p->delta_x[k]);
    gsl_vector_scale(p->mu[k], betak);

    gsl_matrix_memcpy(p->Sigma[k], p->delta_xx[k]);
    gsl_matrix_scale(p->Sigma[k], betak);
    gsl_matrix_view mucv = gsl_matrix_view_array(p->mu[k]->data, p->mu[k]->size, 1);
    gsl_matrix* muc = &(mucv.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, muc, muc, 1.0, p->Sigma[k]);
  }
}

void mcmclib_mixem_online_update(mcmclib_mixem_online* p, gsl_vector* y) {
  (p->n)++;
  mixem_update_s(p, y);
  mixem_update_gamma(p);
}
