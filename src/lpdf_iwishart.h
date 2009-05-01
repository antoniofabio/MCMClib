/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_IWISHART_H__
#define __MCMCLIB_IWISHART_H__

/**\addtogroup distributions
 @{*/
/**\addtogroup multivariate
 @{*/

/**\brief Inverse Wishart distribution */
typedef struct {
  gsl_matrix* Psi; /**< location parameter */
  int m; /**< degrees of freedom */
  double PsiDet; /**< log-determinant of location parameter x m/2*/

  gsl_matrix *PsiX, *X1; /**< workspace data */
} mcmclib_iwishart_lpdf;

/**\brief alloc a Wishart distribution object */
mcmclib_iwishart_lpdf* mcmclib_iwishart_lpdf_alloc(gsl_matrix* Psi, int m);
/**\brief de-alloc a Wishart distribution object */
void mcmclib_iwishart_lpdf_free(mcmclib_iwishart_lpdf* p);

/**\brief compute log-InverseWishart distribution */
double mcmclib_iwishart_lpdf_compute(void* in_p, gsl_vector* x);

/**@}*/
/**@}*/
#endif
