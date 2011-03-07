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

/**@{
\defgroup Base Mean, Variance, Acceptance Rate, Mean Squared Jumping Distance
@{*/

/** Scalar MCMC diagnostics on a 'monitored' vector */
typedef struct mcmclib_monitor_t* mcmclib_monitor_h;

/** alloc a new monitor object */
mcmclib_monitor_h mcmclib_monitor_alloc(const gsl_vector* x);
/** free a previously allocated monitor object */
void mcmclib_monitor_free(mcmclib_monitor_h p);
/** updates monitor statistics using current vector value */
int mcmclib_monitor_update(mcmclib_monitor_h p);

/** get scalar means */
void mcmclib_monitor_get_means(const mcmclib_monitor_h p, gsl_vector* out);
/** get scalar variances */
void mcmclib_monitor_get_vars(const mcmclib_monitor_h p, gsl_vector* out);
/** get scalar acceptance rates */
void mcmclib_monitor_get_ar(const mcmclib_monitor_h p, gsl_vector* out);
/** get scalar MSJD */
void mcmclib_monitor_get_msjd(const mcmclib_monitor_h p, gsl_vector* out);

void mcmclib_monitor_fprintf_means(const mcmclib_monitor_h p, FILE* f);
void mcmclib_monitor_fprintf_vars(const mcmclib_monitor_h p, FILE* f);
void mcmclib_monitor_fprintf_AR(const mcmclib_monitor_h p, FILE* f);
void mcmclib_monitor_fprintf_MSJD(const mcmclib_monitor_h p, FILE* f);
/** formatted print of all the diagnostics on file 'f'*/
void mcmclib_monitor_fprintf_all(const mcmclib_monitor_h p, FILE* f);
/**@}*/
/**@}*/

/**@{
\defgroup ECDF Empirical Cumulative Distrubution Function
@{*/
typedef struct mcmclib_monitor_ecdf_t* mcmclib_monitor_ecdf_h;

/** Alloc a new monitor_ecdf object*/
mcmclib_monitor_ecdf_h mcmclib_monitor_ecdf_alloc(const gsl_matrix* X0);
/** Free a previously allocated monitor_ecdf object*/
void mcmclib_monitor_ecdf_free(mcmclib_monitor_ecdf_h p);
/** Update the ECDF*/
void mcmclib_monitor_ecdf_update(mcmclib_monitor_ecdf_h p, const gsl_vector* y);
/** Get the current ECDF estimate*/
void mcmclib_monitor_ecdf_get(const mcmclib_monitor_ecdf_h p, gsl_vector* Fn);
/**@}*/
/**@}*/


/**@{
\defgroup ACF AutoCorrelation Function
@{*/
typedef struct mcmclib_monitor_acf_t* mcmclib_monitor_acf_h;
/** Alloc a new monitor_acf object*/
mcmclib_monitor_acf_h mcmclib_monitor_acf_alloc(const size_t dim, const size_t lag);
/** Free a previously allocated monitor_acf object*/
void mcmclib_monitor_acf_free(mcmclib_monitor_acf_h m);
/** Update the ACF*/
void mcmclib_monitor_acf_update(mcmclib_monitor_acf_h m, const gsl_vector* x);
/** Get the current ACF estimate*/
void mcmclib_monitor_acf_get(mcmclib_monitor_acf_h m, gsl_matrix* acf);
void mcmclib_iact_from_acf(const gsl_matrix* ACF, gsl_vector* iact);
/**@}*/
/**@}*/

/**@}*/
/**@}*/

#endif
