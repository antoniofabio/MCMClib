%{
#include <monitor.h>
  %}

/** Scalar MCMC diagnostics on a 'monitored' vector */
typedef struct {
  const gsl_vector* x; /**< current value */

  gsl_vector *sum_x, *sum_xsq, *AR, *SJD;
  double n;

  /*internal stuff*/
  gsl_vector *xm, *xvar, *xsq, *ar, *msjd;
  gsl_vector* x_last;
} mcmclib_monitor;

%extend mcmclib_monitor {
  mcmclib_monitor(const gsl_vector* x) {
    return mcmclib_monitor_alloc(x);
  }
  ~mcmclib_monitor() {
    mcmclib_monitor_free($self);
  }
}

int mcmclib_monitor_update(mcmclib_monitor* p);
void mcmclib_monitor_get_means(mcmclib_monitor* p, gsl_vector* out);
void mcmclib_monitor_get_vars(mcmclib_monitor* p, gsl_vector* out);
void mcmclib_monitor_get_ar(mcmclib_monitor* p, gsl_vector* out);
void mcmclib_monitor_get_msjd(mcmclib_monitor* p, gsl_vector* out);

void mcmclib_monitor_update_all(mcmclib_monitor* p);
void mcmclib_monitor_fprintf_means(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_vars(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_AR(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_MSJD(mcmclib_monitor* p, FILE* f);
void mcmclib_monitor_fprintf_all(mcmclib_monitor* p, FILE* f);

/** Empirical CDF function live computation */
typedef struct {
  gsl_vector* Fn; /**< current CDF value */
  gsl_matrix* X0; /**< vertical matrix of points on which the CDF is computed*/
  double n; /**< current observed sample size*/
  gsl_vector* workspace;
} mcmclib_monitor_ecdf;

%extend mcmclib_monitor_ecdf {
  mcmclib_monitor_ecdf(const gsl_matrix* X0) {
    return mcmclib_monitor_ecdf_alloc(X0);
  }
  ~mcmclib_monitor_ecdf() {
    mcmclib_monitor_ecdf_free($self);
  }
}

void mcmclib_monitor_ecdf_update(mcmclib_monitor_ecdf* p, const gsl_vector* y);
