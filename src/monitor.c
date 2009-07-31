/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <gsl/gsl_math.h>
#include "matrix.h"
#include "monitor.h"

mcmclib_monitor* mcmclib_monitor_alloc(const gsl_vector* x) {
  mcmclib_monitor* p = (mcmclib_monitor*) malloc(sizeof(mcmclib_monitor));
  p->x = x;
  int n = x->size;
  p->sum_x = gsl_vector_alloc(n);
  gsl_vector_set_zero(p->sum_x);
  p->sum_xsq = gsl_vector_alloc(n);
  gsl_vector_set_zero(p->sum_xsq);
  p->AR = gsl_vector_alloc(n);
  gsl_vector_set_zero(p->AR);
  p->SJD = gsl_vector_alloc(n);
  gsl_vector_set_zero(p->SJD);

  p->xm = gsl_vector_alloc(n);
  p->xvar = gsl_vector_alloc(n);
  p->xsq = gsl_vector_alloc(n);
  p->msjd = gsl_vector_alloc(n);
  p->ar = gsl_vector_alloc(n);

  p->x_last = gsl_vector_alloc(n);
  gsl_vector_memcpy(p->x_last, x);

  p->n = 0.0;
  return p;
}

void mcmclib_monitor_free(mcmclib_monitor* p) {
  gsl_vector_free(p->x_last);

  gsl_vector_free(p->ar);
  gsl_vector_free(p->xm);
  gsl_vector_free(p->xvar);
  gsl_vector_free(p->xsq);
  gsl_vector_free(p->msjd);

  gsl_vector_free(p->SJD);
  gsl_vector_free(p->AR);
  gsl_vector_free(p->sum_xsq);
  gsl_vector_free(p->sum_x);
  free(p);
}

static void vec_is_null(gsl_vector* x) {
  for(int i=0; i < x->size; i++)
    gsl_vector_set(x, i, gsl_vector_get(x, i) > EQ_TOL ? 1.0 : 0.0);
}

static void vec_sq(gsl_vector* dest, const gsl_vector* x) {
  for(int i=0; i < x->size; i++) {
    double xi = gsl_vector_get(x, i);
    gsl_vector_set(dest, i, xi*xi);
  }
}

int mcmclib_monitor_update(mcmclib_monitor* p) {
  if(!mcmclib_vector_is_finite(p->x))
    GSL_ERROR("non-finite vector value", GSL_EDOM);

  const gsl_vector* x = p->x;

  p->n += 1.0;
  gsl_vector_add(p->sum_x, x);

  vec_sq(p->xsq, x);
  gsl_vector_add(p->sum_xsq, p->xsq);

  gsl_vector_memcpy(p->xsq, p->x_last);
  gsl_vector_memcpy(p->x_last, x);
  gsl_vector_sub(p->xsq, x);
  vec_sq(p->xsq, p->xsq);
  gsl_vector_add(p->SJD, p->xsq);

  vec_is_null(p->xsq);
  gsl_vector_add(p->AR, p->xsq);

  return GSL_SUCCESS;
}

static void update_means(mcmclib_monitor* p) {
  double n1 = 1.0 / p->n;
  gsl_vector_memcpy(p->xm, p->sum_x);
  gsl_vector_scale(p->xm, n1);
}

/** IMPORTANT NOTE: update_means must be called before that func.*/
static void update_variances(mcmclib_monitor* p) {
  double n1 = 1.0 / p->n;
  gsl_vector_memcpy(p->xvar, p->sum_xsq);
  gsl_vector_scale(p->xvar, n1);
  vec_sq(p->xm, p->xm);
  gsl_vector_sub(p->xvar, p->xm);
}

static void update_AR(mcmclib_monitor* p) {
  double n1 = 1.0 / (p->n - 1.0);
  gsl_vector_memcpy(p->ar, p->AR);
  gsl_vector_scale(p->ar, n1);
}

static void update_MSJD(mcmclib_monitor* p) {
  double n1 = 1.0 / (p->n - 1.0);
  gsl_vector_memcpy(p->msjd, p->SJD);
  gsl_vector_scale(p->msjd, n1);
}

void mcmclib_monitor_update_all(mcmclib_monitor* p) {
  update_means(p);
  update_variances(p);
  update_AR(p);
  update_MSJD(p);
}

void mcmclib_monitor_get_means(mcmclib_monitor* p, gsl_vector* out) {
  update_means(p);
  gsl_vector_memcpy(out, p->xm);
}
void mcmclib_monitor_get_vars(mcmclib_monitor* p, gsl_vector* out) {
  update_variances(p);
  gsl_vector_memcpy(out, p->xvar);
}
void mcmclib_monitor_get_ar(mcmclib_monitor* p, gsl_vector* out) {
  update_AR(p);
  gsl_vector_memcpy(out, p->ar);
}
void mcmclib_monitor_get_msjd(mcmclib_monitor* p, gsl_vector* out) {
  update_MSJD(p);
  gsl_vector_memcpy(out, p->msjd);
}

void mcmclib_monitor_fprintf_means(mcmclib_monitor* p, FILE* f) {
  update_means(p);
  gsl_vector_fprintf(f, p->xm, "%f");
}

void mcmclib_monitor_fprintf_vars(mcmclib_monitor* p, FILE* f) {
  update_means(p);
  update_variances(p);
  gsl_vector_fprintf(f, p->xvar, "%f");
}

void mcmclib_monitor_fprintf_AR(mcmclib_monitor* p, FILE* f) {
  update_AR(p);
  gsl_vector_fprintf(f, p->ar, "%f");
}

void mcmclib_monitor_fprintf_MSJD(mcmclib_monitor* p, FILE* f) {
  update_MSJD(p);
  gsl_vector_fprintf(f, p->msjd, "%f");
}

void mcmclib_monitor_fprintf_all(mcmclib_monitor* p, FILE* f) {
  mcmclib_monitor_fprintf_means(p, f);
  mcmclib_monitor_fprintf_vars(p, f);
  mcmclib_monitor_fprintf_AR(p, f);
  mcmclib_monitor_fprintf_MSJD(p, f);
}

mcmclib_monitor_ecdf* mcmclib_monitor_ecdf_alloc(const gsl_matrix* X0) {
  mcmclib_monitor_ecdf* ans = (mcmclib_monitor_ecdf*) malloc(sizeof(mcmclib_monitor_ecdf));
  ans->X0 = gsl_matrix_alloc(X0->size1, X0->size2);
  gsl_matrix_memcpy(ans->X0, X0);
  ans->Fn = gsl_vector_alloc(X0->size1);
  gsl_vector_set_zero(ans->Fn);
  ans->n = 0.0;
  ans->workspace = gsl_vector_alloc(X0->size2);
  return ans;
}
void mcmclib_monitor_ecdf_free(mcmclib_monitor_ecdf* p) {
  gsl_vector_free(p->workspace);
  gsl_vector_free(p->Fn);
  gsl_matrix_free(p->X0);
  free(p);
}
void mcmclib_monitor_ecdf_update(mcmclib_monitor_ecdf* p, const gsl_vector* y) {
  gsl_vector_scale(p->Fn, p->n);
  for(int i=0; i < p->X0->size1; i++) {
    gsl_vector_view rv = gsl_matrix_row(p->X0, i);
    gsl_vector_memcpy(p->workspace, &(rv.vector));
    gsl_vector_sub(p->workspace, y);
    gsl_vector_set(p->Fn, i, gsl_vector_get(p->Fn, i) + (!gsl_vector_ispos(p->workspace)));
  }
  p->n += 1.0;
  gsl_vector_scale(p->Fn, 1.0 / p->n);
}
