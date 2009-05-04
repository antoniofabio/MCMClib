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

/* avoids runtime error, for checking matrix for positive definiteness */
static inline double quiet_sqrt (double x) {
  return (x >= 0) ? sqrt(x) : GSL_NAN;
}

static int
try_cholesky (gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  size_t i,j,k;

  /* Do the first 2 rows explicitly.  It is simple, and faster.  And
   * one can return if the matrix has only 1 or 2 rows.  
   */

  double A_00 = gsl_matrix_get (A, 0, 0);
      
  double L_00 = quiet_sqrt(A_00);
      
  if (A_00 <= 0) {
    return GSL_EDOM;
  }

  gsl_matrix_set (A, 0, 0, L_00);
  
  if (M > 1) {
    double A_10 = gsl_matrix_get (A, 1, 0);
    double A_11 = gsl_matrix_get (A, 1, 1);
          
    double L_10 = A_10 / L_00;
    double diag = A_11 - L_10 * L_10;
    double L_11 = quiet_sqrt(diag);
          
    if (diag <= 0) {
      return GSL_EDOM;
    }

    gsl_matrix_set (A, 1, 0, L_10);        
    gsl_matrix_set (A, 1, 1, L_11);
  }
      
  for (k = 2; k < M; k++) {
    double A_kk = gsl_matrix_get (A, k, k);

    for (i = 0; i < k; i++) {
      double sum = 0;

      double A_ki = gsl_matrix_get (A, k, i);
      double A_ii = gsl_matrix_get (A, i, i);

      gsl_vector_view ci = gsl_matrix_row (A, i);
      gsl_vector_view ck = gsl_matrix_row (A, k);

      if (i > 0) {
	gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
	gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);
                
	gsl_blas_ddot (&di.vector, &dk.vector, &sum);
      }

      A_ki = (A_ki - sum) / A_ii;
      gsl_matrix_set (A, k, i, A_ki);
    } 

    {
      gsl_vector_view ck = gsl_matrix_row (A, k);
      gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);

      double sum = gsl_blas_dnrm2 (&dk.vector);
      double diag = A_kk - sum * sum;

      double L_kk = quiet_sqrt(diag);
            
      if (diag <= 0) {
	return GSL_EDOM;
      }
            
      gsl_matrix_set (A, k, k, L_kk);
    }
  }

  /* Now copy the transposed lower triangle to the upper triangle,
   * the diagonal is common.  
   */
      
  for (i = 1; i < M; i++) {
    for (j = 0; j < i; j++) {
      double A_ij = gsl_matrix_get (A, i, j);
      gsl_matrix_set (A, j, i, A_ij);
    }
  } 

  return GSL_SUCCESS;
}

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
  for(int n= (sigma_chol->size1 - 1); n>=0; n--)
    gsl_vector_set(out, n, gsl_ran_gaussian(r, 1.0));

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
  if(mcmclib_mvnorm_lpdf_chol(p))
    return log(0.0);

  return mcmclib_mvnorm_lpdf_compute_nochol(p, x);
}

int mcmclib_mvnorm_lpdf_chol(mcmclib_mvnorm_lpdf* p) {
  int d = p->mean->size;
  gsl_matrix_view mv = gsl_matrix_view_array(p->vcov, d, d);
  gsl_matrix* vcov = &(mv.matrix);
  gsl_matrix_memcpy(p->rooti, vcov);
  return try_cholesky(p->rooti);
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
