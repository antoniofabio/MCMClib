/*definitions for the target distribution*/
static gsl_vector* mu[K];
static gsl_matrix* Sigma[K];
static mcmclib_mvnorm_lpdf* pi[K];

/*target log-distribution: mixture of two normals*/
double target_logdensity(void* ignore, gsl_vector* x) {
  return log(BETA * exp(mcmclib_mvnorm_lpdf_compute(pi[0], x)) +
	     (1-BETA) * exp(mcmclib_mvnorm_lpdf_compute(pi[1], x)));
}
static void target_distrib_init() {
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
}
static void target_distrib_free() {
  for(int k=0; k<K; k++) {
    gsl_vector_free(mu[k]);
    gsl_matrix_free(Sigma[k]);
    mcmclib_mvnorm_lpdf_free(pi[k]);
  }
}
