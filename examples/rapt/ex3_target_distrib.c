/*definitions for the target distribution*/
static gsl_vector* mu1;
static gsl_vector* mu2;
static gsl_matrix* Sigma1;
static gsl_matrix* Sigma2;
static mcmclib_mvnorm_lpdf* pi1;
static mcmclib_mvnorm_lpdf* pi2;

/*target log-distribution: mixture of two normals*/
double target_logdensity(void* ignore, gsl_vector* x) {
  return log(exp(mcmclib_mvnorm_lpdf_compute(pi1, x)) +
	     exp(mcmclib_mvnorm_lpdf_compute(pi2, x)));
}
static void target_distrib_init() {
  mu1 = gsl_vector_alloc(DIM);
  mu2 = gsl_vector_alloc(DIM);
  Sigma1 = gsl_matrix_alloc(DIM, DIM);
  Sigma2 = gsl_matrix_alloc(DIM, DIM);
  gsl_vector_set_all(mu1, MU0);
  gsl_vector_set_all(mu2, - MU0);

  gsl_matrix_set_identity(Sigma1);
  gsl_matrix_scale(Sigma1, V0);
  gsl_matrix_memcpy(Sigma2, Sigma1);
  /*  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_matrix* tmp1 = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix* tmp2 = gsl_matrix_alloc(DIM, DIM);
  for(int i=0; i<DIM; i++) for(int j=0; j<DIM; j++) {
      gsl_matrix_set(tmp1, i, j, gsl_ran_gaussian(r, 1.0));
      gsl_matrix_set(tmp2, i, j, gsl_ran_gaussian(r, 1.0));
    }
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp1, tmp1, 0.0, Sigma1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp2, tmp2, 0.0, Sigma2);
  gsl_matrix_free(tmp1);
  gsl_matrix_free(tmp2);
  gsl_rng_free(r);*/
  pi1 = mcmclib_mvnorm_lpdf_alloc(mu1, Sigma1->data);
  pi2 = mcmclib_mvnorm_lpdf_alloc(mu2, Sigma2->data);
}
static void target_distrib_free() {
  gsl_vector_free(mu1);
  gsl_vector_free(mu2);
  gsl_matrix_free(Sigma1);
  gsl_matrix_free(Sigma2);
  mcmclib_mvnorm_lpdf_free(pi1);
  mcmclib_mvnorm_lpdf_free(pi2);
}
