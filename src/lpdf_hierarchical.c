/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "lpdf_hierarchical.h"

mcmclib_post_lpdf* mcmclib_post_lpdf_alloc(gsl_vector* x, distrfun_p prior, void* parms,
					   distrfun_p loglik, gsl_vector** childs,
					   void** child_parms, int nchilds) {
  mcmclib_post_lpdf* ans = (mcmclib_post_lpdf*) malloc(sizeof(mcmclib_post_lpdf));
  ans->x = x;
  ans->prior = prior;
  ans->parms = parms;
  ans->loglik = loglik;
  ans->childs = childs;
  ans->child_parms = child_parms;
  ans->nchilds = nchilds;
  ans->workspace = gsl_vector_alloc(x->size);
  return ans;
}

void mcmclib_post_lpdf_free(mcmclib_post_lpdf* p) {
  gsl_vector_free(p->workspace);
  free(p);
}

double mcmclib_post_lpdf_compute(void* data, gsl_vector* x) {
  mcmclib_post_lpdf* d = (mcmclib_post_lpdf*) data;
  double ans = 0.0;
  /*store old value*/
  gsl_vector_memcpy(d->workspace, d->x);
  /*set new value*/
  gsl_vector_memcpy(d->x, x);

  ans += d->prior(d->parms, d->x);
  if(!isfinite(ans)) {
    /*restore old value*/
    gsl_vector_memcpy(d->x, d->workspace);
    return(ans);
  }
  for(int i=0; i < d->nchilds; i++) {
    ans += d->loglik(d->child_parms[i], d->childs[i]);
    if(!isfinite(ans))
      break;
  }

  /*restore old value*/
  gsl_vector_memcpy(d->x, d->workspace);

  return ans;
}
