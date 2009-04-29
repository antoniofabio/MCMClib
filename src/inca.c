/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "inca.h"

mcmclib_inca* mcmclib_inca_alloc(mcmclib_amh* amh, gsl_vector** x, int M) {
  mcmclib_inca* p = (mcmclib_inca*) malloc(sizeof(mcmclib_inca));
  p->amh = amh;
  p->x = x;
  p->M = M;
  return p;
}

void mcmclib_inca_free(mcmclib_inca* p) {
  free(p);
}

int mcmclib_inca_update(mcmclib_inca* p) {
  mcmclib_amh* amh = p->amh;
  mcmclib_mh* mh = amh->mh;
  gsl_vector** x = p->x;
  for(int m=0; m < p->M; m++) {
    gsl_vector_memcpy(mh->x, x[m]);
    mcmclib_amh_update(amh);
    gsl_vector_memcpy(x[m], mh->x);
  }
  return(mh->last_accepted);
}
