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
#include <memory.h>
#include "mixnorm.h"

mcmclib_mixnorm_lpdf* mcmclib_mixnorm_lpdf_alloc(gsl_vector* w, mcmclib_mvnorm_lpdf** pis){
  mcmclib_mixnorm_lpdf* a = (mcmclib_mixnorm_lpdf*) malloc(sizeof(mcmclib_mixnorm_lpdf));
  a->pis = (mcmclib_mvnorm_lpdf**) malloc(w->size * sizeof(mcmclib_mvnorm_lpdf*));
  memcpy(a->pis, pis, w->size * sizeof(mcmclib_mvnorm_lpdf*));
  for(size_t i=0; i<w->size; i++)
    a->pis[i] = pis[i];
  a->w = gsl_vector_alloc(w->size);
  gsl_vector_memcpy(a->w, w);
  return a;
}

void mcmclib_mixnorm_lpdf_free(mcmclib_mixnorm_lpdf* p) {
  free(p->pis);
  gsl_vector_free(p->w);
  free(p);
}

double mcmclib_mixnorm_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mixnorm_lpdf* p = (mcmclib_mixnorm_lpdf*) in_p;
  double ans = 0.0;
  for(size_t i = 0; i < p->w->size; i++)
    ans += exp(mcmclib_mvnorm_lpdf_compute(p->pis[i], x)) * gsl_vector_get(p->w, i);
  return log(ans);
}
