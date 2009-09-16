%module mcmclib

%include "carrays.i"

%scheme %{(load-extension "libguilemcmclib.so" "scm_init_swig_mcmclib_module")%}

%{
#include <lpdf_iwishart.h>
#include <mh.h>
#include <monitor.h>
#include <gauss_rw.h>
#include <gauss_mrw.h>

static void guile_mcmclib_err_handler(const char * reason,
				      const char * file,
				      int line,
				      int gsl_errno) {
  static char msg[2048];
  sprintf(msg, "%s:%d", file, line);
  scm_error_scm(scm_misc_error_key,
		scm_from_locale_string(msg),
		scm_from_locale_string(reason), SCM_BOOL_F, SCM_BOOL_F);
}
%}

%init %{
  gsl_set_error_handler(guile_mcmclib_err_handler);
%}

typedef struct {
  gsl_rng* r; /**< rng*/
  distrfun_p logdistr; /**< target log-density fun*/
  void* logdistr_data; /**< target log-density data*/
  mcmclib_mh_q* q; /**< proposal kernel object*/
  gsl_vector* x; /**< current chain value*/
  gsl_vector* x_old; /**< old chain value*/
  double logdistr_old; /**< log-distrib. of old chain value */
  int last_accepted; /**< flag: last move has been accepted?*/
} mcmclib_mh;
int mcmclib_mh_update(mcmclib_mh* p);

typedef double (*distrfun_p)(void*, gsl_vector*);

mcmclib_mh* mcmclib_gauss_rw_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* data,
				   gsl_vector* start_x, double step_size);
void mcmclib_gauss_rw_free(mcmclib_mh*);
mcmclib_mh* mcmclib_gauss_mrw_alloc(gsl_rng* r,
				    distrfun_p distrfun,
				    void* logdistr_data,
				    gsl_vector* x,
				    const gsl_matrix* sigma_zero);
void mcmclib_gauss_mrw_free(mcmclib_mh*);

mcmclib_iwishart_lpdf* mcmclib_iwishart_lpdf_alloc(gsl_matrix* Psi, int m);
void mcmclib_iwishart_lpdf_free(mcmclib_iwishart_lpdf* p);
%callback("%s_cb");
double mcmclib_iwishart_lpdf_compute(void* p, gsl_vector* x);
%nocallback;

%{
  void* test_distrfun_alloc(double what) {
    double* ans = malloc(sizeof(double));
    ans[0] = what;
    return ans;
  }

  double test_distrfun(void* p, gsl_vector* x) {
    double* pwhat = (double*) p;
    double x0 = gsl_vector_get(x, 0);
    if((x0 >= 0.0) && (x0 <= pwhat[0]))
      return log(1.0);
    else
      return log(0.0);
  }
%}
void* test_distrfun_alloc(double what);
%callback("%s_cb");
double test_distrfun(void* p, gsl_vector* x);
%nocallback;

%{
  static double mcmclib_guile_lpdf(void* p, gsl_vector* x) {
    SCM sx = SWIG_NewPointerObj(x, SWIGTYPE_p_gsl_vector, 0);
    SCM ans = scm_call_1((SCM) p, sx);
    return scm_to_double(ans);
  }
%}
%inline %{
  void* guile_to_voidptr(SCM p) {
    return (void*) p;
  }
%}

%callback("%s_cb");
double mcmclib_guile_lpdf(void* p, gsl_vector* x);
%nocallback;

%include "monitor.i"
%include "distrfuns.i"
%include "amh.i"
%include "mixem.i"
