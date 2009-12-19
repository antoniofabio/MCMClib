/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "mh_q.h"

mcmclib_mh_q* mcmclib_mh_q_alloc(gsl_rng* r,
				 sampler_fun_t rq,
				 proposal_distr_fun_t dq,
				 void* gamma, free_fun_t free_gamma_fun) {
  mcmclib_mh_q* a = (mcmclib_mh_q*) malloc(sizeof(mcmclib_mh_q));
  a->r = r;
  a->rq = rq;
  a->dq = dq;
  a->gamma = gamma;
  a->free_gamma_fun = free_gamma_fun;
  return a;
}

void mcmclib_mh_q_free(mcmclib_mh_q* p) {
  if(p->free_gamma_fun)
    p->free_gamma_fun(p->gamma);
  free(p);
}

void mcmclib_mh_q_sample(mcmclib_mh_q* p, gsl_vector* x) {
  p->rq(p, x);
}

double mcmclib_mh_q_logd(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y) {
  return p->dq(p, x, y);
}

double mcmclib_mh_q_ratio_offset(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y) {
  return mcmclib_mh_q_logd(p, y, x) - mcmclib_mh_q_logd(p, x, y);
}
