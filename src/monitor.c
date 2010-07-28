/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <assert.h>
#include <gsl/gsl_math.h>
#include "vector_queue.h"
#include "matrix.h"
#include "monitor.h"

mcmclib_monitor* mcmclib_monitor_alloc(const gsl_vector* x) {
  mcmclib_monitor* p = (mcmclib_monitor*) malloc(sizeof(mcmclib_monitor));
  p->x = x;
  const size_t n = x->size;
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
  for(size_t i=0; i < x->size; i++)
    gsl_vector_set(x, i, gsl_vector_get(x, i) > EQ_TOL ? 1.0 : 0.0);
}

static void vec_sq(gsl_vector* dest, const gsl_vector* x) {
  for(size_t i=0; i < x->size; i++) {
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
  for(size_t i=0; i < p->X0->size1; i++) {
    gsl_vector_view rv = gsl_matrix_row(p->X0, i);
    gsl_vector_memcpy(p->workspace, &(rv.vector));
    gsl_vector_sub(p->workspace, y);
    gsl_vector_set(p->Fn, i, gsl_vector_get(p->Fn, i) + (!gsl_vector_ispos(p->workspace)));
  }
  p->n += 1.0;
  gsl_vector_scale(p->Fn, 1.0 / p->n);
}

struct mcmclib_monitor_acf_t {
  mcmclib_vector_queue* q;
  gsl_vector_uint* x_n; /*num. obs. so far, by lag*/
  gsl_matrix* x_sum;
  gsl_matrix* X_prod;
};

static size_t monitor_acf_maxlag(const mcmclib_monitor_acf_h m) {
  return m->X_prod->size1 - 1;
}

static size_t monitor_acf_dim(const mcmclib_monitor_acf_h m) {
  return m->x_sum->size2;
}

mcmclib_monitor_acf_h mcmclib_monitor_acf_alloc(const size_t dim, const size_t maxlag) {
  mcmclib_monitor_acf_h m = (mcmclib_monitor_acf_h) malloc(sizeof(struct mcmclib_monitor_acf_t));
  m->q = mcmclib_vector_queue_alloc(dim, maxlag+1);
  m->x_n = gsl_vector_uint_alloc(maxlag+1);
  gsl_vector_uint_set_zero(m->x_n);
  m->x_sum = gsl_matrix_alloc(maxlag+1, dim);
  gsl_matrix_set_zero(m->x_sum);
  m->X_prod = gsl_matrix_alloc(maxlag+1, dim);
  gsl_matrix_set_zero(m->X_prod);
  return m;
}

void mcmclib_monitor_acf_free(mcmclib_monitor_acf_h m) {
  if(!m) return;
  mcmclib_vector_queue_free(m->q);
  gsl_vector_uint_free(m->x_n);
  gsl_matrix_free(m->x_sum);
  gsl_matrix_free(m->X_prod);
  free(m);
}

#define DEBUG_PRINT_VEC(x) fprintf(stderr, "(%f, ...)\n", gsl_vector_get(x, 0))

#define MIN(a, b) ((a) < (b) ? a : b)

void mcmclib_monitor_acf_update(mcmclib_monitor_acf_h m, const gsl_vector* x) {
  const size_t dim = monitor_acf_dim(m);
  assert(x->size == dim);
  fprintf(stderr, "acf_update: appending vector ");
  DEBUG_PRINT_VEC(x);

  const size_t maxlag = monitor_acf_maxlag(m);

  mcmclib_vector_queue_append(m->q, x);

  gsl_vector* xy = gsl_vector_alloc(dim);
  gsl_vector* y = gsl_vector_alloc(dim);
  const size_t qsize = mcmclib_vector_queue_size(m->q);
  for(size_t l = 0; l < qsize; l++) {
    fprintf(stderr, "maxlag = %zd, lag %zd\n", maxlag, l);
    mcmclib_vector_queue_get(m->q, l, y);

    gsl_vector_uint_set(m->x_n, l, gsl_vector_uint_get(m->x_n, l)+1);

    gsl_vector_view x_sum_l_v = gsl_matrix_row(m->x_sum, l);
    gsl_vector_add(&(x_sum_l_v.vector), y);
    fprintf(stderr, "acf_update: lag %zd, partial lagged sum ", l);
    DEBUG_PRINT_VEC(&(x_sum_l_v.vector));

    gsl_vector_memcpy(xy, x);
    gsl_vector_mul(xy, y);
    gsl_vector_view xl_v = gsl_matrix_row(m->X_prod, l);
    gsl_vector_add(&xl_v.vector, xy);
    fprintf(stderr, "acf_update: lag %zd, partial crossed product sum ", l);
    DEBUG_PRINT_VEC(&xl_v.vector);
  }
  gsl_vector_free(y);
  gsl_vector_free(xy);
}

#define VECTOR_MAP(x, op) for (size_t i = 0; i < x->size; i++) {	\
    double xi = gsl_vector_get(x, i);					\
    op;									\
    gsl_vector_set(x, i, xi);						\
}

void mcmclib_monitor_acf_get(mcmclib_monitor_acf_h m, gsl_matrix* acf) {
  const size_t dim = monitor_acf_dim(m);
  const size_t qsize = mcmclib_vector_queue_size(m->q);
  gsl_vector* x_sq = gsl_vector_alloc(dim);

  for(size_t l = 0; l < qsize; l++) {
    gsl_vector_view x_cov_v = gsl_matrix_row(acf, l);
    gsl_vector* x_cov = &(x_cov_v.vector);
    gsl_vector_const_view xx_l_v = gsl_matrix_const_row(m->X_prod, l);
    gsl_vector_memcpy(x_cov, &(xx_l_v.vector));
    gsl_vector_const_view x0_sum_v = gsl_matrix_const_row(m->x_sum, l);
    gsl_vector_memcpy(x_sq, &(x0_sum_v.vector));
    VECTOR_MAP(x_sq, xi*=xi);
    gsl_vector_sub(x_cov, x_sq);
  }

  gsl_vector_view var_v = gsl_matrix_row(acf, 0);
  gsl_vector* var = &(var_v.vector);
  assert(!gsl_vector_isneg(var));
  for(size_t l = 1; l < qsize; l++) {
    gsl_vector_view cv_l_v = gsl_matrix_row(acf, l);
    gsl_vector* cv_l = &(cv_l_v.vector);
    VECTOR_MAP(cv_l, xi /= gsl_vector_get(var, i));
  }

  gsl_vector_free(x_sq);

  gsl_matrix_scale(acf, 1.0 / ((double) gsl_vector_uint_get(m->x_n, 0)));
}

void mcmclib_iact_from_acf(const gsl_matrix* ACF, gsl_vector* iact) {
  const size_t dim = ACF->size2;
  const size_t L = ACF->size1;
  assert(iact->size == dim);
  for(size_t d = 0; d < dim; d++) {
    double a = 1.0;
    for(size_t l = 1; l <= L; l++) {
      a += (1.0 - (double) l / ((double) L)) * gsl_matrix_get(ACF, l, d);
    }
    gsl_vector_set(iact, d, a);
  }
}
