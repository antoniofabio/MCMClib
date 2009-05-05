/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <math.h>
#include <gsl/gsl_math.h>
#include "mvnorm.h"

void mcmclib_mvnorm(const gsl_rng* r,
		    const gsl_matrix* sigma,
		    gsl_vector* out) {
  int d = sigma->size1;

  gsl_matrix* sigma_chol = gsl_matrix_alloc(d, d);
  gsl_matrix_memcpy(sigma_chol, sigma);
  gsl_linalg_cholesky_decomp(sigma_chol);

  mcmclib_mvnorm_chol(r, sigma_chol, out);

  gsl_matrix_free(sigma_chol);
}

void mcmclib_mvnorm_chol(const gsl_rng* r,
			 const gsl_matrix* sigma_chol,
			 gsl_vector* out) {
  /*generate d iid values*/
  for(int n=0; n < sigma_chol->size1; n++)
    gsl_vector_set(out, sigma_chol->size1 - n - 1, gsl_ran_gaussian(r, 1.0));

  /*rotate them according to the cholesky 'square root' of sigma*/
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, sigma_chol, out);
}

mcmclib_mvnorm_lpdf* mcmclib_mvnorm_lpdf_alloc(gsl_vector* mean, double* vcov) {
  int d = mean->size;
  mcmclib_mvnorm_lpdf* ans = (mcmclib_mvnorm_lpdf*) malloc(sizeof(mcmclib_mvnorm_lpdf));
  ans->mean = mean;
  ans->vcov = vcov;
  ans->rooti = gsl_matrix_alloc(d, d);
  ans->x_mu = gsl_vector_alloc(d);
  ans->mahal = gsl_vector_alloc(d);
  return ans;
}

void mcmclib_mvnorm_lpdf_free(mcmclib_mvnorm_lpdf* p) {
  gsl_matrix_free(p->rooti);
  gsl_vector_free(p->x_mu);
  gsl_vector_free(p->mahal);
  free(p);
}

double mcmclib_mvnorm_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mvnorm_lpdf* p = (mcmclib_mvnorm_lpdf*) in_p;

  /*compute cholesky decomposition of var/cov matrix*/
  if(mcmclib_mvnorm_lpdf_chol(p) != GSL_SUCCESS)
    return log(0.0);

  return mcmclib_mvnorm_lpdf_compute_nochol(p, x);
}

int mcmclib_mvnorm_lpdf_chol(mcmclib_mvnorm_lpdf* p) {
  int d = p->mean->size;
  gsl_matrix_view mv = gsl_matrix_view_array(p->vcov, d, d);
  gsl_matrix* vcov = &(mv.matrix);
  gsl_matrix_memcpy(p->rooti, vcov);
  gsl_error_handler_t *hnd = gsl_set_error_handler_off();
  int ans = gsl_linalg_cholesky_decomp(p->rooti);
  gsl_set_error_handler(hnd);
  return ans;
}

double mcmclib_mvnorm_lpdf_compute_nochol(mcmclib_mvnorm_lpdf* p, gsl_vector* x) {
  int d = x->size;

  /*compute mahlanobis distance between 'x' and 'mu'*/
  gsl_vector* x_mu = p->x_mu;
  gsl_vector_memcpy(x_mu, x);
  gsl_vector_sub(x_mu, p->mean);
  gsl_vector_memcpy(p->mahal, x_mu);
  gsl_linalg_cholesky_svx(p->rooti, p->mahal);

  /*compute log-density as:
    -0.5 * (mahaldist + log(2*pi)*d + logdet) */
  double ans = 0.0;
  gsl_blas_ddot(p->mahal, x_mu, &ans);
  ans += log(2.0 * M_PI) * ((double) d);
  for(int i=0; i<d; i++)
    ans += log(gsl_matrix_get(p->rooti, i, i)) * 2.0;
  ans *= -0.5;

  return ans;
}

void mcmclib_mvnorm_lpdf_inverse(mcmclib_mvnorm_lpdf* p) {
  int d = p->mean->size;
  gsl_matrix* xx = p->rooti;
  gsl_matrix_view vcov_v = gsl_matrix_view_array(p->vcov, d, d);
  gsl_matrix_memcpy(xx, &(vcov_v.matrix));
  gsl_permutation* perm = gsl_permutation_alloc(d);
  int signum = 0;
  gsl_linalg_LU_decomp(xx, perm, &signum);
  p->determinant = log(gsl_linalg_LU_det(xx, signum));
  gsl_matrix* tmp = gsl_matrix_alloc(d, d);
  gsl_matrix_memcpy(tmp, xx);
  gsl_linalg_LU_invert(tmp, perm, xx);

  gsl_matrix_free(tmp);
  gsl_permutation_free(perm);
}

double mcmclib_mvnorm_lpdf_compute_noinv(mcmclib_mvnorm_lpdf* p, gsl_vector* x) {
  int d = x->size;

  /*compute mahlanobis distance between 'x' and 'mu'*/
  gsl_vector* x_mu = p->x_mu;
  gsl_vector_memcpy(x_mu, x);
  gsl_vector_sub(x_mu, p->mean); /*x - mu*/
  /*mahal = inv(Sigma) (x - mu)*/
  gsl_blas_dgemv(CblasNoTrans, 1.0, p->rooti, p->x_mu, 0.0, p->mahal);
  double ans = 0.0;
  /*ans = t(mahal) * (x - mu)*/
  gsl_blas_ddot(p->mahal, x_mu, &ans);

  /*compute log-density as:
    -0.5 * (mahaldist + log(2*pi)*d + logdet) */
  ans += log(2.0 * M_PI) * ((double) d) + p->determinant;
  ans *= -0.5;

  return ans;
}
