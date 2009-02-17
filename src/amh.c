#include "amh.h"

mcmclib_amh* mcmclib_amh_alloc(mcmclib_mh* mh,
			       mcmclib_amh_update_gamma_p update_gamma,
			       void* update_gamma_data) {
  mcmclib_amh* p = (mcmclib_amh*) malloc(sizeof(mcmclib_amh));
  p->mh = mh;
  p->update_gamma = update_gamma;
  p->update_gamma_data = update_gamma_data;
  p->n = 0;
  return p;
}

void mcmclib_amh_free(mcmclib_amh* p) {
  free(p);
}

int mcmclib_amh_update(mcmclib_amh* p) {
  mcmclib_mh_update(p->mh);
  p->n++;
  p->update_gamma(p->update_gamma_data, p->mh->x);
  return(p->mh->last_accepted);
}
