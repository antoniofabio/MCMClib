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
/**\addtogroup distributions
 @{*/

/**\addtogroup mcar
 @{*/

/**\brief MCAR(\tilde B, \Gamma) distribution */
typedef struct {
  int p; /**< dimension */
  int n; /**< number of points */

  gsl_vector* alpha1; /**< Givens angles for P1 */
  gsl_vector* alpha2; /**< Givens angles for P2 */
  gsl_vector* sigma; /**< B_tilde singular values (p) */
  gsl_matrix* B_tilde; /**< variance par. matrix (p x p) */
	gsl_matrix* Gamma; /**< 'variance of variance' par. matrix (p x p) */

  gsl_matrix* M; /**< adiancency matrix (n x n)*/
  gsl_vector* m; /**< adiancency weights (n)*/

  /*internal stuff*/
  gsl_vector* mu;
  gsl_matrix* vcov;
  mcmclib_mvnorm_lpdf* mvnorm; /**< normal density object */
} mcmclib_mcar_tilde_lpdf;

/** Alloc extra data for an mcar_tilde distribution
    @param p dimension
    @param n spatial locations
    @param M adiancency matrix (n x n)
*/
mcmclib_mcar_tilde_lpdf* mcmclib_mcar_tilde_lpdf_alloc(int p, int n, gsl_matrix* M);

void mcmclib_mcar_tilde_lpdf_set_alpha(mcmclib_mcar_tilde_lpdf* p,
				       gsl_vector* alpha1, gsl_vector* alpha2);
void mcmclib_mcar_tilde_lpdf_set_sigma(mcmclib_mcar_tilde_lpdf* p,
				       gsl_vector* sigma);

/** Free extra data for an mcar_tilde distribution
    @param p
*/
void mcmclib_mcar_tilde_lpdf_free(mcmclib_mcar_tilde_lpdf* p);

/** mcar_tilde log-distribution
    @param in_p extra data, allocated via \ref mcmclib_mcar_tilde_lpdf_alloc
    @return log-pdf
*/
double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x);

/** update current vcov matrix value */
void mcmclib_mcar_tilde_lpdf_update_vcov(mcmclib_mcar_tilde_lpdf* p);

/**@}*/
/**@}*/
#endif
