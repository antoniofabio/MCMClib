/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_MCAR_TILDE_H__
#define __MCMCLIB_MCAR_TILDE_H__

#include "mvnorm.h"
#include "givens.h"
/**\addtogroup distributions
 @{*/

/**\addtogroup mcar
 @{*/

/**\brief MCAR(\tilde B, \Gamma) distribution */
typedef struct {
  int p; /**< dimension */
  int n; /**< number of points */

  gsl_matrix* B_tilde; /**< variance par. matrix (p x p) */
  gsl_vector* alpha12sigma; /**< Givens angles and sing. values repr. of
			       B_tilde */
  gsl_matrix* Gamma; /**< 'variance of variance' par. matrix (p x p) */
  gsl_vector *alphasigmag; /**< Givens angles and eigenv. repr. of Gamma */

  gsl_matrix* M; /**< adiancency matrix (n x n)*/
  gsl_vector* m; /**< adiancency weights (n)*/

  /*internal stuff*/
  gsl_vector* mu;
  gsl_matrix* vcov;
  mcmclib_mvnorm_lpdf* mvnorm; /**< normal density object */

  /*workspace memory*/
  gsl_matrix *Lambda_ij, *Gammai, *Block; /**< used in fun 'vcov_blockij' */
} mcmclib_mcar_tilde_lpdf;

/** Alloc extra data for an mcar_tilde distribution
    @param p dimension
    @param M adiancency matrix (n x n)
*/
mcmclib_mcar_tilde_lpdf* mcmclib_mcar_tilde_lpdf_alloc(int p, gsl_matrix* M);

/** Free extra data for an mcar_tilde distribution
    @param p
*/
void mcmclib_mcar_tilde_lpdf_free(mcmclib_mcar_tilde_lpdf* p);

/** mcar_tilde log-distribution
    @param in_p extra data, allocated via \ref mcmclib_mcar_tilde_lpdf_alloc
    @return log-pdf
*/
double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x);

/** update current vcov matrix value \internal */
void mcmclib_mcar_tilde_lpdf_update_B_tilde(mcmclib_mcar_tilde_lpdf* p);
/** update current inverse vcov matrix value \internal */
void mcmclib_mcar_tilde_lpdf_update_blocks(mcmclib_mcar_tilde_lpdf* p);
/** update current vcov matrix value \internal */
int mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p);
/** update current Gamma matrix value \internal */
void mcmclib_mcar_tilde_lpdf_update_Gamma(mcmclib_mcar_tilde_lpdf* p);

/**@}*/
/**@}*/
#endif
