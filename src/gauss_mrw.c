/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <gsl/gsl_math.h>
#include "gauss_mrw.h"
#include "mvnorm.h"

double mcmclib_gauss_mrw_qd(mcmclib_mh_q* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

void mcmclib_gauss_mrw_sample(mcmclib_mh_q* q, gsl_vector* x) {
  gsl_matrix* Sigma = (gsl_matrix*) q->gamma;
  gsl_vector* x_old = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(x_old, x);
  mcmclib_mvnorm(q->r, Sigma, x);
  gsl_vector_add(x, x_old);
  gsl_vector_free(x_old);
}

void mcmclib_gauss_mrw_free(void* p) {
  gsl_matrix_free((gsl_matrix*) p);
}

mcmclib_mh* mcmclib_gauss_mrw_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* data,
				    gsl_vector* start_x, const gsl_matrix* Sigma) {
  gsl_matrix* Sigma1 = gsl_matrix_alloc(Sigma->size1, Sigma->size2);
  gsl_matrix_memcpy(Sigma1, Sigma);
  return mcmclib_mh_alloc(r, logdistr, data,
			  mcmclib_mh_q_alloc(r, mcmclib_gauss_mrw_sample,
					     mcmclib_gauss_mrw_qd,
					     Sigma1,
					     mcmclib_gauss_mrw_free),
			  start_x);
}
