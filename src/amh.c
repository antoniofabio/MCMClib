#include "amh.h"

mcmclib_amh* mcmclib_amh_alloc(mcmclib_mh* mh, void* suff,
			       mcmclib_amh_update_gamma_p update_gamma) {
  mcmclib_amh* p = (mcmclib_amh*) malloc(sizeof(mcmclib_amh));
  p->mh = mh;
  p->suff = suff;
  p->update_gamma = update_gamma;
  p->n = 0;
  return p;
}

void mcmclib_amh_free(mcmclib_amh* p) {
  free(p);
}

int mcmclib_amh_update(mcmclib_amh* p) {
  mcmclib_mh_update(p->mh);
  p->n++;
  p->update_gamma(p, p->mh->x);
  return(p->mh->last_accepted);
}
