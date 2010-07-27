/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "amh.h"

mcmclib_amh* mcmclib_amh_alloc(mcmclib_mh* mh, const size_t n0,
			       void* suff,
			       free_fun_t free_suff,
			       mcmclib_amh_update_suff_t update_suff,
			       mcmclib_amh_update_gamma_t update_gamma) {
  mcmclib_amh* p = (mcmclib_amh*) malloc(sizeof(mcmclib_amh));
  p->mh = mh;
  p->suff = suff;
  p->update_suff = update_suff;
  p->update_gamma = update_gamma;
  p->n = 0;
  p->n0 = n0;
  p->free_suff = free_suff;
  return p;
}

void mcmclib_amh_reset(mcmclib_amh* p) {
  mcmclib_mh_reset(p->mh);
  p->n = 0;
}

void mcmclib_amh_free(mcmclib_amh* p) {
  if(p->free_suff)
    p->free_suff(p->suff);
  mcmclib_mh_free(p->mh);
  free(p);
}

int mcmclib_amh_update(mcmclib_amh* p) {
  mcmclib_mh_update(p->mh);
  p->n++;
  if(p->update_suff) {
    p->update_suff(p->suff, p->mh->x);
  }
  if(p->update_gamma && (p->n >= p->n0)) {
    p->update_gamma(p->suff, p->mh->q->gamma);
  }
  return(p->mh->last_accepted);
}
