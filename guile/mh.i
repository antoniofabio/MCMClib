%{
#include <mh.h>
#include <gauss_rw.h>
#include <gauss_mrw.h>

  double mcmclib_guile_lpdf(void* p, const gsl_vector* x) {
    SCM sx = SWIG_NewPointerObj(x, SWIGTYPE_p_gsl_vector, 0);
    SCM ans = scm_call_1((SCM) p, sx);
    return scm_to_double(ans);
  }
%}

typedef double (*distrfun_p) (void* data, const gsl_vector* x);

%typemap(in) (distrfun_p f, void* data) {
  scm_permanent_object($input); /*FIXME. Maybe solve it by proper use of '$owner'*/
  $1 = mcmclib_guile_lpdf;
  $2 = (void*) $input;
}

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

mcmclib_mh* mcmclib_mh_alloc(gsl_rng* r, distrfun_p f, void* data,
			     mcmclib_mh_q* q, gsl_vector* x);
void mcmclib_mh_free(mcmclib_mh* p);

int mcmclib_mh_update(mcmclib_mh* p);

typedef double (*distrfun_p)(void*, gsl_vector*);

mcmclib_mh* mcmclib_gauss_rw_alloc(gsl_rng* r,
				   distrfun_p f, void* data,
				   gsl_vector* start_x, double step_size);
mcmclib_mh* mcmclib_gauss_mrw_alloc(gsl_rng* r,
				    distrfun_p f, void* data,
				    gsl_vector* x,
				    const gsl_matrix* sigma_zero);
