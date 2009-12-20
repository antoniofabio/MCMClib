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

void mcmclib_gauss_rw_sample(mcmclib_mh_q* q, gsl_vector* x) {
  int d = x->size;
  double step_size = *((double*)q->gamma);
  while(d--) {
    gsl_vector_set(x, d,
		   gsl_vector_get(x, d) +
		   gsl_ran_gaussian(q->r, step_size));
  }
}
double mcmclib_gauss_rw_qd(mcmclib_mh_q* q, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

mcmclib_mh* mcmclib_gauss_rw_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* data,
				   gsl_vector* start_x, double step_size) {
  return mcmclib_mh_alloc(r, logdistr, data,
			  mcmclib_mh_q_alloc(r,
					     mcmclib_gauss_rw_sample,
					     mcmclib_gauss_rw_qd,
					     &step_size,
					     NULL),
			  start_x);
}
