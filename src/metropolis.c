#include "metropolis.h"

int mcmclib_metropolis_symmetric_step(const gsl_rng* r, gsl_vector* old, gsl_vector* x, distrfun_p logdistr, void* data) {
  double loglik_old, loglik_new, lik_ratio;

  loglik_old = logdistr(data, old);
  if(!isfinite(loglik_old))
    return 1;

  loglik_new = logdistr(data, x);
  if(isfinite(loglik_new) && (loglik_new >= loglik_old))
    return 1;

  lik_ratio = exp(loglik_new - loglik_old);
  if(isfinite(loglik_new) && isfinite(lik_ratio) && (gsl_rng_uniform(r) <= lik_ratio))
    return 1;

  gsl_vector_memcpy(x, old);
  return 0;
}

int mcmclib_metropolis_generic_step(const gsl_rng* r, gsl_vector* old,
				    gsl_vector* x, distrfun_p logdistr, void* data,
				    proposal_distr_t q, void* q_data) {
  double loglik_old, loglik_new;
  double q_old_new, q_new_old;
  double mh_ratio;

  loglik_old = logdistr(data, old);
  if(!isfinite(loglik_old))
    return 1;
  q_old_new = q(q_data, old, x);
  if(!isfinite(q_old_new))
    return 1;

  loglik_new = logdistr(data, x);
  if(!isfinite(loglik_new)) {
    gsl_vector_memcpy(x, old);
    return 0;
  }
  q_new_old = q(q_data, x, old);
  if(!isfinite(q_new_old)) {
    gsl_vector_memcpy(x, old);
    return 0;
  }

  mh_ratio = loglik_new - loglik_old +  q_new_old - q_old_new;
  if((mh_ratio > 0) || (gsl_rng_uniform(r) <= exp(mh_ratio)))
    return 1;

  gsl_vector_memcpy(x, old);
  return 0;
}
