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
#include "gauss_rw.h"

mcmclib_gauss_rw_gamma* mcmclib_gauss_rw_gamma_alloc(gsl_rng* r, double step_size) {
  mcmclib_gauss_rw_gamma* a = (mcmclib_gauss_rw_gamma*) malloc(sizeof(mcmclib_gauss_rw_gamma));
  a->r = r;
  a->step_size = step_size;
  return a;
}
void mcmclib_gauss_rw_gamma_free(mcmclib_gauss_rw_gamma* p) {
  free(p);
}
void mcmclib_gauss_rw_sample(void* in_q, gsl_vector* x) {
  mcmclib_mh_q* q = (mcmclib_mh_q*) in_q;
  mcmclib_gauss_rw_gamma* p = (mcmclib_gauss_rw_gamma*) q->sampler_data;
  int d = x->size;
  while(d--) {
    gsl_vector_set(x, d,
		   gsl_vector_get(x, d) + gsl_ran_gaussian(p->r, p->step_size));
  }
}
double mcmclib_gauss_rw_qd(void* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

mcmclib_mh_q* mcmclib_gauss_rw_q_alloc(gsl_rng* r, double step_size) {
  mcmclib_gauss_rw_gamma* gamma = mcmclib_gauss_rw_gamma_alloc(r, step_size);
  return mcmclib_mh_q_alloc(r, mcmclib_gauss_rw_sample, gamma,
			    mcmclib_gauss_rw_qd, gamma,
			    gamma);
}
void mcmclib_gauss_rw_q_free(mcmclib_mh_q* p) {
  mcmclib_gauss_rw_gamma_free(p->gamma);
  mcmclib_mh_q_free(p);
}

mcmclib_mh* mcmclib_gauss_rw_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* data,
				   gsl_vector* start_x, double step_size) {
  return mcmclib_mh_alloc(r, logdistr, data, mcmclib_gauss_rw_q_alloc(r, step_size),
		   start_x);
}

void mcmclib_gauss_rw_free(mcmclib_mh* p) {
  mcmclib_gauss_rw_q_free(p->q);
  mcmclib_mh_free(p);
}
