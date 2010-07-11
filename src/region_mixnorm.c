/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "mixnorm.h"
#include "region_mixnorm.h"

size_t mcmclib_region_mixnorm_compute(const gsl_vector* x, void* in_p) {
  mcmclib_mixnorm_lpdf* p = (mcmclib_mixnorm_lpdf*) in_p;
  double pik;
  int ans = 0;
  double pimax = mcmclib_mvnorm_lpdf_compute(p->pis[0], x);
  for(size_t k=1; k < p->w->size; k++) {
    pik = mcmclib_mvnorm_lpdf_compute(p->pis[k], x);
    if(pik > pimax) {
      pimax = pik;
      ans = k;
    }
  }
  return ans;
}
