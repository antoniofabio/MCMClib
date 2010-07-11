/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_WISHART_H__
#define __MCMCLIB_WISHART_H__

/**\addtogroup distributions
 @{*/
/**\addtogroup multivariate
 @{*/

/**\brief Wishart distribution */
typedef struct {
  gsl_matrix* V; /**< (inverse of) location parameter */
  unsigned int p; /**< degrees of freedom */

  gsl_matrix *Vx, *X1; /**< workspace data */
} mcmclib_wishart_lpdf;

/**\brief alloc a Wishart distribution object */
mcmclib_wishart_lpdf* mcmclib_wishart_lpdf_alloc(const gsl_matrix* V, const unsigned int p);
/**\brief de-alloc a Wishart distribution object */
void mcmclib_wishart_lpdf_free(mcmclib_wishart_lpdf* p);

/**\brief compute log-Wishart distribution */
double mcmclib_wishart_lpdf_compute(void* in_p, const gsl_vector* x);

/**@}*/
/**@}*/
#endif
