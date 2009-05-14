/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "matrix.h"
#include "lm.h"

mcmclib_lm* mcmclib_lm_alloc() {
  mcmclib_lm* ans;
  return ans;
}
void mcmclib_lm_free(mcmclib_lm* p) {
  free(p);
}

/*ported from the MCMCpack R package C++ sources*/
int NormNormregress_beta_draw(const gsl_matrix* XX, const gsl_vector* XY,
			      const gsl_vector* b0, const gsl_matrix* B0,
			      double sigma2, gsl_rng* rng, gsl_vector* out) {
  // this function gets the cross-product matrix X'X and the matrix X'Y
  // to minimize the amount of computation within the function
  assert(XX->size1 == XX->size2);
  const unsigned int k = XX->size2;
  const double sig2_inv = 1.0 / sigma2;
  gsl_matrix* sig_beta = gsl_matrix_alloc(k, k);
  gsl_matrix_memcpy(sig_beta, XX);
  gsl_matrix_scale(sig_beta, sig2_inv);
  gsl_matrix_add(sig_beta, B0);
  mcmclib_cholesky_inverse(sig_beta);
  gsl_matrix* C = gsl_matrix_alloc(k, k);
  gsl_matrix_memcpy(C, sig_beta);
  gsl_linalg_cholesky_decomp(C);
  gsl_vector* betahat = gsl_vector_alloc(k);
  gsl_vector_memcpy(betahat, XY);
  //const Matrix<> betahat = sig_beta * gaxpy(B0, b0, XpY*sig2_inv);
  gsl_blas_dgemv(CblasNoTrans, 1.0, B0, b0, sig2_inv, betahat);
  gsl_blas_dgemv(CblasNoTrans, 1.0, sig_beta, betahat, 0.0, betahat);
  mcmclib_mvnorm_chol(rng, C, out);
  gsl_vector_add(out, betahat);
  gsl_vector_free(betahat);
  gsl_matrix_free(C);
  gsl_matrix_free(sig_beta);
  return GSL_SUCCESS;
}

/*ported from the MCMCpack R package C++ sources*/
// linear regression with Gaussian errors sigma2 draw 
// (inverse-Gamma  prior)
// regression model is y = X * beta + epsilon,  epsilon ~ N(0,sigma2)
// c0/2 is the prior shape parameter for sigma2
// d0/2 is the prior scale parameter for sigma2 
int NormIGregress_sigma2_draw (const gsl_matrix* X, const gsl_vector* Y,
			       const gsl_vector* beta, double c0, double d0,
			       gsl_rng* rng, gsl_vector* out) {
  int k = X->size2;
  gsl_vector* e = gsl_vector_alloc(k);
  gsl_vector_memcpy(e, Y);
  gsl_blas_dgemv(CblasNoTrans, -1.0, X, beta, 1.0, e);
  double SSE;
  gsl_blas_ddot(e, e, &SSE);
  gsl_vector_free(e);
  int n = X->size1;
  const double c_post = (c0 + n) * 0.5;
  const double d_post = (d0 + SSE) * 0.5;

  gsl_vector_set(out, 0, 1.0 / gsl_ran_gamma(rng, c_post, d_post));
  return GSL_SUCCESS;
}
