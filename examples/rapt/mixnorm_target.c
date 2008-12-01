#include <mixnorm.h>
#include "mixnorm_target.h"

/*state space dimension*/
const int DIM = 5;
/*absolute local mean value*/
const double MU0 = 1.5;
/*local variances multipliers*/
const double V0[] = {4.0, 1.0};
/*local pairwise correlations*/
const double RHO[] = {-0.1, -0.1};
/*mixture proportion of component 1*/
const double BETA = 0.2;

double beta[K];
gsl_vector_view beta_vv;
gsl_vector* mu[K];
gsl_matrix* Sigma[K];
mcmclib_mvnorm_lpdf* pi[K];

mcmclib_mixnorm_lpdf* mixnorm_target_alloc() {
  beta[0] = BETA;
  beta[1] = 1-BETA;
  for(int k=0; k<K; k++) {
    mu[k] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(mu[k], pow(-1.0, k + 1) * MU0);
    Sigma[k] = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix_set_identity(Sigma[k]);
    gsl_matrix_scale(Sigma[k], V0[k]);
    for(int i=1; i<DIM; i++) for(int j=0; j<i; j++){
	gsl_matrix_set(Sigma[k], i, j, RHO[k]);
	gsl_matrix_set(Sigma[k], j, i, RHO[k]);
    }
    pi[k] = mcmclib_mvnorm_lpdf_alloc(mu[k], Sigma[k]->data);
  }
  beta_vv = gsl_vector_view_array(beta, K);
  return mcmclib_mixnorm_lpdf_alloc(&(beta_vv.vector), pi);
}

void mixnorm_target_free(mcmclib_mixnorm_lpdf* p) {
  for(int k=0; k<K; k++) {
    gsl_vector_free(mu[k]);
    gsl_matrix_free(Sigma[k]);
    mcmclib_mvnorm_lpdf_free(pi[k]);
  }
  mcmclib_mixnorm_lpdf_free(p);
}
