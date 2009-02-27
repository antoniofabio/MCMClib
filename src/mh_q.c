#include "mh.h"

typedef struct {
  gsl_rng* r; /**< RNG used by the sampler function*/
  samplerfun_p sampler_fun; /**< proposal sampler fun*/
  void* sampler_data; /**< proposal sampler data*/
  proposal_distr_t qd_fun; /**< proposal density fun*/
  void* qd_data; /**< proposal density data*/
  void* gamma; /**< misc kernel parameters data*/
} mcmclib_mh_q;

mcmclib_mh_q* mcmclib_mh_q_alloc(gsl_rng* r,
				 samplerfun_p sampler_fun, void* sampler_data,
				 proposal_distr_t qd_fun, void* qd_data,
				 void* gamma) {
  mcmclib_mh_q* a = (mcmclib_mh_q*) malloc(sizeof(mcmclib_mh_q*));
  a->r = r;
  a->sampler_fun = sampler_fun;
  a->sampler_data = sampler_data;
  a->qd_fun = qd_fun;
  a->qd_data = qd_data;
  a->gamma = gamma;
  return a;
}

void mcmclib_mh_q_free(mcmclib_mh_q* p) {
  free(p);
}

void mcmclib_mh_q_sample(mcmclib_mh_q* p, gsl_vector* x) {
  p->sampler_fun(p->sampler_data, x);
}

double mcmclib_mh_q_logd(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y) {
  return p->qd_fun(p->qd_data, x, y);
}

double mcmclib_mh_ratio_offset(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y) {
  return mcmclib_mh_q_logd(p, x, y) - mcmclib_mh_q_logd(p, y, x);
}
