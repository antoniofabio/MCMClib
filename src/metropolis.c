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
  double loglik_old, loglik_new, loglik_diff;
  double q_old_new, q_new_old, q_diff;
  double mh_ratio;

  loglik_old = logdistr(data, old);
  if(!isfinite(loglik_old))
    return 1;
  loglik_new = logdistr(data, x);

  mh_ratio = q(q_data, old, x) - q(q_data, x, old) + loglik_new + loglik_old;
  if(isfinite(mh_ratio) && (gsl_rng_uniform(r) <= mh_ratio))
    return 1;

  gsl_vector_memcpy(x, old);
  return 0;
}
