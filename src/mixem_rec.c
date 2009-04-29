/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include "mixem_rec.h"
#include "mvnorm.h"
#include "mixnorm.h"
#include "vector_stats.h"

mcmclib_mixem_rec* mcmclib_mixem_rec_alloc(gsl_vector** mu,
					   gsl_matrix** Sigma,
					   gsl_vector* beta){
  mcmclib_mixem_rec* a = (mcmclib_mixem_rec*) malloc(sizeof(mcmclib_mixem_rec));
  a->mu = mu;
  a->Sigma = Sigma;
  a->beta = beta;

  a->n = 0;
  int K = beta->size;
  a->pi_k = (mcmclib_mvnorm_lpdf**) malloc(sizeof(mcmclib_mvnorm_lpdf*) * K);
  a->X_sq_sum = (gsl_matrix**) malloc(sizeof(gsl_matrix*) * K);
  a->X_sum = (gsl_vector**) malloc(sizeof(gsl_vector*) * K);
  for(int k=0; k<K; k++) {
    a->pi_k[k] = mcmclib_mvnorm_lpdf_alloc(mu[k], Sigma[k]->data);
    a->X_sq_sum[k] = gsl_matrix_alloc(mu[k]->size, mu[k]->size);
    gsl_matrix_set_all(a->X_sq_sum[k], 0.0);
    a->X_sum[k] = gsl_vector_alloc(mu[k]->size);
    gsl_vector_set_all(a->X_sum[k], 0.0);
  }
  a->beta_sum = gsl_vector_alloc(beta->size);
  gsl_vector_set_all(a->beta_sum, 0.0);
  a->beta_i = gsl_vector_alloc(beta->size);

  return a;
}

void mcmclib_mixem_rec_free(mcmclib_mixem_rec* p) {
  for(int k=0; k < p->beta->size; k++) {
    mcmclib_mvnorm_lpdf_free(p->pi_k[k]);
    gsl_matrix_free(p->X_sq_sum[k]);
    gsl_vector_free(p->X_sum[k]);
  }
  gsl_vector_free(p->beta_sum);
  gsl_vector_free(p->beta_i);
  free(p->X_sq_sum);
  free(p->X_sum);
  free(p->pi_k);
  free(p);
}

void mcmclib_mixem_rec_add(mcmclib_mixem_rec* p, gsl_vector* y) {
  int K = p->beta->size; /*# mixture components*/
  (p->n)++;

  /*store posterior class probs. of y in beta_i*/
  double pi_sum = 0.0;
  for(int k = 0; k<K; k++) {
    double pik = exp(mcmclib_mvnorm_lpdf_compute(p->pi_k[k], y)) * gsl_vector_get(p->beta, k);
    gsl_vector_set(p->beta_i, k, pik);
    pi_sum += pik;
  }
  gsl_vector_scale(p->beta_i, 1.0 / pi_sum);

  /*update weighted sum of p(yi|theta), yi, yi*t(yi) for each class*/
  gsl_vector_add(p->beta_sum, p->beta_i);
  gsl_matrix_view myv = gsl_matrix_view_array(y->data, y->size, 1);
  gsl_matrix* my = &(myv.matrix);
  for(int k=0; k<K; k++) {
    double wik = gsl_vector_get(p->beta_i, k);
    gsl_blas_daxpy(wik, y, p->X_sum[k]);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,
		   wik, my, my, 1.0, p->X_sq_sum[k]);
  }
}

/*update theta estimates basing on current accumulated values*/
void mcmclib_mixem_rec_update(mcmclib_mixem_rec* p) {
  int K = p->beta->size;
  if(p->n < 2)
    return;
  for(int k=0; k<K; k++) {
    double wik = 1.0 / gsl_vector_get(p->beta_sum, k);

    gsl_vector_memcpy(p->mu[k], p->X_sum[k]);
    gsl_vector_scale(p->mu[k], wik);

    gsl_matrix_view mmuv = gsl_matrix_view_array(p->mu[k]->data, p->mu[k]->size, 1);
    gsl_matrix* mmu = &(mmuv.matrix);

    gsl_matrix_memcpy(p->Sigma[k], p->X_sq_sum[k]);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, mmu, mmu, wik, p->Sigma[k]);
  }
  gsl_vector_memcpy(p->beta, p->beta_sum);
  gsl_vector_scale(p->beta, 1.0 / (double) p->n);
}
