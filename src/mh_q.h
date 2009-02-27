#ifndef __MH_Q_H__
#define __MH_Q_H__

/**\addtogroup metropolis_samplers
@{
\defgroup MH Metropolis-Hastings sampling
*/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

/** pointer to a sampler function */
typedef void (*samplerfun_p) (void* data, gsl_vector* x);

/**\brief Metropolis-Hastings proposal kernel*/
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
				 void* gamma);
void mcmclib_mh_q_free(mcmclib_mh_q* p);

/**\brief sample a new proposal point*/
void mcmclib_mh_q_sample(mcmclib_mh_q* p, gsl_vector* x);
/**\brief wrapper around the proposal kernel log-density function*/
double mcmclib_mh_q_logd(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y);
/**\brief compute the M-H ratio offset using the kernel q function*/
double mcmclib_mh_ratio_offset(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y);

/**@}*/
/**@}*/

#endif
