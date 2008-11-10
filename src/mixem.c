#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "mixem.h"
#include "mvnorm.h"
#include "mixnorm.h"
#include "vector_stats.h"

/**Fitting a gaussian mixture distribution by the EM algorithm
*/
void mcmclib_mixem_fit(gsl_matrix* X, int K,
		       gsl_vector** mu, gsl_matrix** Sigma,
		       gsl_matrix* P, gsl_vector* w, int NITER) {
  mcmclib_mvnorm_lpdf** pi_k = (mcmclib_mvnorm_lpdf**) malloc(K * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++) {
    pi_k[k] = mcmclib_mvnorm_lpdf_alloc(mu[k], Sigma[k]->data);
    mcmclib_mvnorm_lpdf_chol(pi_k[k]);
  }
  int N = X->size1;
  int d = X->size2;
  gsl_vector* wn = w;
  gsl_vector* tmp_X_n = gsl_vector_alloc(d);

  for(int ITER=0; ITER < NITER; ITER++) {
    /*get weight prediction for each point*/
    gsl_vector_set_all(wn, 0.0);
    for(int n=0; n<N; n++) {
      gsl_vector_view rv = gsl_matrix_row(X, n);
      gsl_vector* r = &(rv.vector);
      double rowsum = 0.0;
      for(int k=0; k<K; k++) {
	double p_nk = exp(mcmclib_mvnorm_lpdf_compute_nochol(pi_k[k], r));
	gsl_matrix_set(P, n, k, p_nk);
	rowsum += p_nk;
      }
      gsl_vector_view P_n = gsl_matrix_row(P, n);
      gsl_vector_scale(&(P_n.vector), 1.0 / rowsum);
      gsl_vector_add(w, &(P_n.vector));
    }
    
    /*update means and variances estimates*/
    for(int k=0; k<K; k++) {
      gsl_vector_set_all(mu[k], 0.0);
      gsl_matrix_set_all(Sigma[k], 0.0);
      gsl_vector_view P_kv = gsl_matrix_column(P, k);
      gsl_vector* P_k = &(P_kv.vector);
      for(int n=0; n<N; n++) {
	double P_nk = gsl_vector_get(P_k, n);
	gsl_vector_view X_nv = gsl_matrix_row(X, n);
	gsl_vector* X_n = &(X_nv.vector);
	gsl_vector_memcpy(tmp_X_n, X_n);
	gsl_vector_scale(tmp_X_n, P_nk);
	gsl_vector_add(mu[k], tmp_X_n);
	gsl_vector_memcpy(tmp_X_n, X_n);
	gsl_matrix_view X_n_colv = gsl_matrix_view_array(tmp_X_n->data, d, 1);
	gsl_matrix* X_n_col = &(X_n_colv.matrix);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, P_nk,
		       X_n_col, X_n_col, 1.0, Sigma[k]);
      }
      gsl_vector_scale(mu[k], 1.0 / gsl_vector_get(w, k));
      gsl_matrix_view mu_colv = gsl_matrix_view_array(mu[k]->data, d, 1);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0,
		     &(mu_colv.matrix), &(mu_colv.matrix),
		     1.0 / gsl_vector_get(w, k), Sigma[k]);
    }
    gsl_vector_scale(w, 1.0 / (double) N);

  }

  gsl_vector_free(tmp_X_n);
  for(int k=0; k<K; k++)
    mcmclib_mvnorm_lpdf_free(pi_k[k]);
  free(pi_k);
}
