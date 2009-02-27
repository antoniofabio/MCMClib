#include "mh.h"

mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r,
			     distrfun_p logdistr, void* logdistr_data,
			     mcmclib_mh_q* q, gsl_vector* x) {
  mcmclib_mh* p = (mcmclib_mh*) malloc(sizeof(mcmclib_mh));
  p->r = r;
  p->logdistr = logdistr;
  p->logdistr_data = logdistr_data;
  p->q = q;
  p->x = x;
  p->x_old = gsl_vector_alloc(x->size);
  p->last_accepted = 0;
  return p;
}

void mcmclib_mh_free(mcmclib_mh* p) {
  gsl_vector_free(p->x_old);
  free(p);
}

int mcmclib_mh_update(mcmclib_mh* p) {
  gsl_vector_memcpy(p->x_old, p->x);
  mcmclib_mh_q_sample(p->q, p->x);
  p->last_accepted = mcmclib_mh_generic_step(p->r, p->x_old, p->x,
					     p->logdistr, p->logdistr_data,
					     p->q);
  return(p->last_accepted);
}

int mcmclib_mh_generic_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x,
			    distrfun_p logdistr, void* data, mcmclib_mh_q* q) {
  double logdistr_old, logdistr_new;
  double mh_offset, mh_ratio;

  logdistr_old = logdistr(data, old);
  if(!isfinite(logdistr_old))
    return 1;
  mh_offset = mcmclib_mh_q_ratio_offset(q, old, x);
  if(!isfinite(mh_offset)) {
    if(mh_offset < 0) {
      gsl_vector_memcpy(x, old);
      return 0;
    } else
      return 1;
  }

  logdistr_new = logdistr(data, x);
  if(!isfinite(logdistr_new)) {
    gsl_vector_memcpy(x, old);
    return 0;
  }

  mh_ratio = logdistr_new - logdistr_old +  mh_offset;
  if((mh_ratio > 0) || (gsl_rng_uniform(r) <= exp(mh_ratio)))
    return 1;

  gsl_vector_memcpy(x, old);
  return 0;
}
