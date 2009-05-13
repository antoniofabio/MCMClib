/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MONITOR_H__
#define __MONITOR_H__

/**\addtogroup misc
@{
\defgroup Monitor General MCMC live diagnostics
@{*/

#include <gsl/gsl_vector.h>

#define EQ_TOL 1e-10 /**< equality check tolerance*/

/** Scalar MCMC diagnostics on a 'monitored' vector */
typedef struct {
  const gsl_vector* x; /**< current value */

  gsl_vector *sum_x, *sum_xsq, *AR, *SJD;
  double n;

  /*internal stuff*/
  gsl_vector *xm, *xvar, *xsq, *ar, *msjd;
  gsl_vector* x_last;
} mcmclib_monitor;

/** alloc a new monitor object */
mcmclib_monitor* mcmclib_monitor_alloc(const gsl_vector* x);
/** free a previously allocated monitor object */
void mcmclib_monitor_free(mcmclib_monitor* p);
/** updates monitor statistics using current vector value */
int mcmclib_monitor_update(mcmclib_monitor* p);

/** get scalar means */
void mcmclib_monitor_get_means(mcmclib_monitor* p, gsl_vector* out);
/** get scalar variances */
void mcmclib_monitor_get_vars(mcmclib_monitor* p, gsl_vector* out);
/** get scalar acceptance rates */
void mcmclib_monitor_get_ar(mcmclib_monitor* p, gsl_vector* out);
/** get scalar MSJD */
void mcmclib_monitor_get_msjd(mcmclib_monitor* p, gsl_vector* out);

void mcmclib_monitor_update_all(mcmclib_monitor* p);
void mcmclib_monitor_fprintf_means(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_vars(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_AR(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_MSJD(mcmclib_monitor* p, FILE* f);
/** formatted print of all the diagnostics on file 'f'*/
void mcmclib_monitor_fprintf_all(mcmclib_monitor* p, FILE* f);

/**@}*/
/**@}*/

#endif
