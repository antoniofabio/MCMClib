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

mcmclib_amh* mcmclib_amh_alloc(mcmclib_mh* mh, void* suff,
			       mcmclib_amh_update_gamma_p update_gamma,
			       free_fun_t free_suff) {
  mcmclib_amh* p = (mcmclib_amh*) malloc(sizeof(mcmclib_amh));
  p->mh = mh;
  p->suff = suff;
  p->update_gamma = update_gamma;
  p->n = 0;
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
  p->update_gamma(p);
  return(p->mh->last_accepted);
}
