#ifndef __RAPTOR_H__
#define __RAPTOR_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "mixem_online.h"
#include "mixnorm.h"
#include "rapt.h"

/**\addtogroup adaptive
@{*/
/**\defgroup RAPTOR raptor
\brief RAPT based on On-Line EM fitting of a Gaussian mixture
@{*/

/** \brief RAPTOR sampler gamma values */
typedef struct {
  gsl_vector* beta_hat; /**< current mixture weights estimates*/
  gsl_vector** mu_hat; /**< current mixture means estimates*/
  gsl_matrix** Sigma_hat; /**< current mixture variances estimates*/

  mcmclib_mvnorm_lpdf** pik_hat; /**< single mixture components densities*/
  mcmclib_mixnorm_lpdf* pi_hat; /**< mixture density*/
} mcmclib_raptor_gamma;

/** \brief RAPTOR sufficient data */
typedef struct {
  mcmclib_mixem_online* em; /**< online-EM mixture fitter*/
} mcmclib_raptor_suff;

/** alloc a new RAPTOR sampler suff. stats. object
@param t0 burn-in length before starting adaptation
@returns a new raptor_suff object
*/
mcmclib_raptor_suff* mcmclib_raptor_suff_alloc(mcmclib_raptor_gamma* g, int t0);

/** free raptor_suff data*/
void mcmclib_raptor_suff_free(mcmclib_raptor_suff* p);

/** Update current value of an OLEM-RAPT chain*/
int mcmclib_raptor_suff_update(mcmclib_raptor_suff* p);

/** update local and global RAPT proposals covariance matrices

basing on current mixture parameters estimates*/
void mcmclib_raptor_update_proposals(mcmclib_raptor_suff* p);

/**@}*/
/**@}*/
#endif
