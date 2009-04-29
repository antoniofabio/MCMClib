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

mcmclib_gauss_mrw_gamma* mcmclib_gauss_mrw_gamma_alloc(gsl_rng* r, const gsl_matrix* Sigma) {
  mcmclib_gauss_mrw_gamma* a = (mcmclib_gauss_mrw_gamma*) malloc(sizeof(mcmclib_gauss_mrw_gamma));
  a->r = r;
  a->Sigma = gsl_matrix_alloc(Sigma->size1, Sigma->size1);
  gsl_matrix_memcpy(a->Sigma, Sigma);
  return a;
}
void mcmclib_gauss_mrw_gamma_free(mcmclib_gauss_mrw_gamma* p) {
  gsl_matrix_free(p->Sigma);
  free(p);
}
void mcmclib_gauss_mrw_sample(void* in_p, gsl_vector* x) {
  mcmclib_mh_q* q = (mcmclib_mh_q*) in_p;
  mcmclib_gauss_mrw_gamma* p = (mcmclib_gauss_mrw_gamma*) q->gamma;
  gsl_vector* x_old = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(x_old, x);
  mcmclib_mvnorm(p->r, p->Sigma, x);
  gsl_vector_add(x, x_old);
  gsl_vector_free(x_old);
}
double mcmclib_gauss_mrw_qd(void* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

mcmclib_mh_q* mcmclib_gauss_mrw_q_alloc(gsl_rng* r, const gsl_matrix* Sigma) {
  mcmclib_gauss_mrw_gamma* gamma = mcmclib_gauss_mrw_gamma_alloc(r, Sigma);
  return mcmclib_mh_q_alloc(r, mcmclib_gauss_mrw_sample, gamma,
			    mcmclib_gauss_mrw_qd, gamma,
			    gamma);
}

void mcmclib_gauss_mrw_q_free(mcmclib_mh_q* p) {
  mcmclib_gauss_mrw_gamma_free(p->gamma);
  mcmclib_mh_q_free(p);
}

mcmclib_mh* mcmclib_gauss_mrw_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* data,
				    gsl_vector* start_x, const gsl_matrix* Sigma) {
  return mcmclib_mh_alloc(r, logdistr, data, mcmclib_gauss_mrw_q_alloc(r, Sigma),
		   start_x);
}

void mcmclib_gauss_mrw_free(mcmclib_mh* p) {
  mcmclib_gauss_mrw_q_free(p->q);
  mcmclib_mh_free(p);
}
