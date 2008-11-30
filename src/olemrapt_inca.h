#ifndef __OLEMRAPT_INCA_H__
#define __OLEMRAPT_INCA_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "olem_rapt.h"

/**\addtogroup adaptive
@{*/
/**\defgroup OLEMRAPT_INCA olemrapt_inca
\brief OLEM-RAPT Inter-Chain Adaptation
@{*/

/** OLEMRAPT-INCA data */
typedef struct {
  int M; /**< number of parallel samplers*/
  mcmclib_olemrapt** ss; /**< array of olemrapt samplers*/

  mcmclib_mixolem_suff* s; /**< pooled sufficient statistic*/
  mcmclib_mixolem_suff* workspace; /**< workspace memory*/
  mcmclib_mixolem_suff* gamma_hat; /**< mixture parameters estimates*/
} mcmclib_olemrapt_inca;

/** alloc a new OLEM-RAPT INCA sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x array of current chains values
@param t0 burn-in length before starting adaptation
@param Sigma_zero initial global proposal covariance matrix
@param beta_hat starting weights estimates
@param mu_hat starting means estimates
@param Sigma_hat starting variances estimates
@return a new olemrapt_inca object
*/
mcmclib_olemrapt_inca* mcmclib_olemrapt_inca_alloc(gsl_rng* r,
						   distrfun_p logdistr, void* logdistr_data,
						   gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
						   gsl_vector* beta_hat,
						   gsl_vector** mu_hat,
						   gsl_matrix** Sigma_hat);

/** free OLEM-RAPT INCA data*/
void mcmclib_olemrapt_inca_free(mcmclib_olemrapt_inca* p);

/** Update current values of OLEM-RAPT INCA chains*/
int mcmclib_olemrapt_inca_update(mcmclib_olemrapt_inca* p);

/** update local and global RAPT proposals covariance matrices

basing on current mixture parameters pooled estimates*/
void mcmclib_olemrapt_inca_update_proposals(mcmclib_olemrapt_inca* p);

/**@}*/
/**@}*/
#endif
