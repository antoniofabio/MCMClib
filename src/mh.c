#include "mh.h"

mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r,
			     distrfun_p logdistr,
			     void* logdistr_data,
			     gsl_vector* x,
			     proposal_distr_t qd_fun,
			     void* qd_data,
			     samplerfun_p q_sampler,
			     void* q_sampler_data) {
  mcmclib_mh* p = (mcmclib_mh*) malloc(sizeof(mcmclib_mh));
  p->r = r;
  p->logdistr = logdistr;
  p->logdistr_data = logdistr_data;
  p->x = x;
  p->x_old = gsl_vector_alloc(x->size);
  p->qd_fun = qd_fun;
  p->qd_data = qd_data;
  p->q_sampler = q_sampler;
  p->q_sampler_data = q_sampler_data;
  p->last_accepted = 0;
  return p;
}

void mcmclib_mh_free(mcmclib_mh* p) {
  gsl_vector_free(p->x_old);
  free(p);
}

int mcmclib_mh_update(mcmclib_mh* p) {
  gsl_vector_memcpy(p->x_old, p->x);
  p->q_sampler(p->q_sampler_data, p->x);
  p->last_accepted = mcmclib_metropolis_generic_step(p->r, p->x_old, p->x,
						     p->logdistr, p->logdistr_data,
						     p->qd_fun, p->qd_data);
  return(p->last_accepted);
}
