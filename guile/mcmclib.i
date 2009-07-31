%module mcmclib

%include "carrays.i"

%scheme %{(load-extension "libguilemcmclib.so" "scm_init_swig_mcmclib_module")%}

%{
#include <lpdf_iwishart.h>
#include <mh.h>
#include <monitor.h>
#include <gauss_mrw.h>
%}

int mcmclib_mh_update(mcmclib_mh* p);

typedef double (*distrfun_p)(void*, gsl_vector*);

%newobject mcmclib_gauss_mrw_alloc;
mcmclib_mh* mcmclib_gauss_mrw_alloc(gsl_rng* r,
				    distrfun_p distrfun,
				    void* logdistr_data,
				    gsl_vector* x,
				    const gsl_matrix* sigma_zero);
%delobject mcmclib_gauss_mrw_free;
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
