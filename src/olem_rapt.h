#ifndef __OLEMRAPT_H__
#define __OLEMRAPT_H__

/**\file
\brief OnLine EM-based RAPT
*/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "mixem_online.h"
#include "mixnorm.h"
#include "rapt.h"

/** OLEM-RAPT data */
typedef struct {
  mcmclib_rapt* rapt;

  gsl_vector* beta_hat;
  gsl_vector** mu_hat;
  gsl_matrix** Sigma_hat;

  mcmclib_mvnorm_lpdf** pik_hat;
  mcmclib_mixnorm_lpdf* pi_hat;

  mcmclib_mixem_online* em;
} mcmclib_olemrapt;

/** alloc a new OLEM-RAPT sampler object
@param r RNG state
@param logdistr pointer to a log-likelihood function
@param logdistr_data extra data to be passed to the distribution function
@param x current chain value
@param t0 burn-in length before starting adaptation
@param Sigma_zero initial global proposal covariance matrix
@param beta_hat starting weights estimates
@param mu_hat starting means estimates
@param Sigma_hat starting variances estimates
*/
mcmclib_olemrapt* mcmclib_olemrapt_alloc(gsl_rng* r,
					 distrfun_p logdistr, void* logdistr_data,
					 gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
					 gsl_vector* beta_hat,
					 gsl_vector** mu_hat,
					 gsl_matrix** Sigma_hat);

/** free OLEM-RAPT data*/
void mcmclib_olemrapt_free(mcmclib_olemrapt* p);

/** Update current value of an OLEM-RAPT chain*/
int mcmclib_olemrapt_update(mcmclib_olemrapt* p);

/** update local and global proposals covariance matrices

basing on current mixture parameters estimates*/
void mcmclib_olemrapt_update_proposals(mcmclib_olemrapt* p);

#endif
