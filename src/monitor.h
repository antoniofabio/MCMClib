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

/** update all the monitored stats*/
void mcmclib_monitor_update_all(mcmclib_monitor* p);
void mcmclib_monitor_fprintf_means(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_vars(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_AR(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_MSJD(mcmclib_monitor* p, FILE* f);
/** formatted print of all the diagnostics on file 'f'*/
void mcmclib_monitor_fprintf_all(mcmclib_monitor* p, FILE* f);

typedef struct mcmclib_monitor_ecdf_t* mcmclib_monitor_ecdf_h;

/** Empirical CDF function live computation */
mcmclib_monitor_ecdf_h mcmclib_monitor_ecdf_alloc(const gsl_matrix* X0);
void mcmclib_monitor_ecdf_free(mcmclib_monitor_ecdf_h p);
void mcmclib_monitor_ecdf_update(mcmclib_monitor_ecdf_h p, const gsl_vector* y);

typedef struct mcmclib_monitor_acf_t* mcmclib_monitor_acf_h;
/** Autocorrelation function live computation */
mcmclib_monitor_acf_h mcmclib_monitor_acf_alloc(const size_t dim, const size_t lag);
void mcmclib_monitor_acf_update(mcmclib_monitor_acf_h m, const gsl_vector* x);
void mcmclib_monitor_acf_free(mcmclib_monitor_acf_h m);
void mcmclib_monitor_acf_get(mcmclib_monitor_acf_h m, gsl_matrix* acf);
void mcmclib_iact_from_acf(const gsl_matrix* ACF, gsl_vector* iact);

/**@}*/
/**@}*/

#endif
