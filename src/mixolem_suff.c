/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "mixolem_suff.h"

mcmclib_mixolem_suff* mcmclib_mixolem_suff_makeref(gsl_vector* delta,
						   gsl_vector** delta_x,
						   gsl_matrix** delta_xx) {
  mcmclib_mixolem_suff* a = (mcmclib_mixolem_suff*) malloc(sizeof(mcmclib_mixolem_suff));
  a->delta = delta;
  a->delta_x = delta_x;
  a->delta_xx = delta_xx;
  return a;
}

mcmclib_mixolem_suff* mcmclib_mixolem_suff_alloc(size_t K, size_t dim) {
  gsl_vector* delta = gsl_vector_alloc(K);
  gsl_vector_set_all(delta, 0.0);
  gsl_vector** delta_x = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  gsl_matrix** delta_xx = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(size_t k=0; k < K; k++) {
    delta_x[k] = gsl_vector_alloc(dim);
    gsl_vector_set_all(delta_x[k], 0.0);
    delta_xx[k] = gsl_matrix_alloc(dim, dim);
    gsl_matrix_set_all(delta_xx[k], 0.0);
  }
  return mcmclib_mixolem_suff_makeref(delta, delta_x, delta_xx);
}

void mcmclib_mixolem_suff_free(mcmclib_mixolem_suff* p) {
  for(size_t k=0; k < p->delta->size; k++) {
    gsl_vector_free(p->delta_x[k]);
    gsl_matrix_free(p->delta_xx[k]);
  }
  free(p->delta_x);
  free(p->delta_xx);
  gsl_vector_free(p->delta);
  free(p);
}

void mcmclib_mixolem_suff_memcpy(mcmclib_mixolem_suff* dest, mcmclib_mixolem_suff* src){
  gsl_vector_memcpy(dest->delta, src->delta);
  for(size_t k=0; k < dest->delta->size; k++) {
    gsl_vector_memcpy(dest->delta_x[k], src->delta_x[k]);
    gsl_matrix_memcpy(dest->delta_xx[k], src->delta_xx[k]);
  }
}

void mcmclib_mixolem_suff_add(mcmclib_mixolem_suff* dest, mcmclib_mixolem_suff* src) {
  gsl_vector_add(dest->delta, src->delta);
  for(size_t k=0; k < dest->delta->size; k++) {
    gsl_vector_add(dest->delta_x[k], src->delta_x[k]);
    gsl_matrix_add(dest->delta_xx[k], src->delta_xx[k]);
  }
}

void mcmclib_mixolem_suff_scale(mcmclib_mixolem_suff* p, double alpha) {
  gsl_vector_scale(p->delta, alpha);
  for(size_t k=0; k < p->delta->size; k++) {
    gsl_vector_scale(p->delta_x[k], alpha);
    gsl_matrix_scale(p->delta_xx[k], alpha);
  }
}
