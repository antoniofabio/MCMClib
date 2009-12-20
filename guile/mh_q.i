%{
#include <mh_q.h>
  %}

struct mcmclib_mh_q_t;

%{
  double mcmclib_mh_q_guile_dq(struct mcmclib_mh_q_t* data,
			       gsl_vector* x, gsl_vector* y) {
    SCM gamma = (SCM) data->gamma;
    SCM dq = scm_cadr(gamma);
    SCM sx = scm_gc_protect_object(SWIG_NewPointerObj(x, SWIGTYPE_p_gsl_vector, 0));
    SCM sy = scm_gc_protect_object(SWIG_NewPointerObj(y, SWIGTYPE_p_gsl_vector, 0));
    double ans = scm_to_double(scm_call_2(dq, sx, sy));
    scm_gc_unprotect_object(sx);
    scm_gc_unprotect_object(sy);
    return ans;
  }
  void mcmclib_mh_q_guile_rq(struct mcmclib_mh_q_t* data,
			       gsl_vector* x) {
    SCM gamma = (SCM) data->gamma;
    SCM rq = scm_car(gamma);
    SCM sx = scm_gc_protect_object(SWIG_NewPointerObj(x, SWIGTYPE_p_gsl_vector, 0));
    scm_call_1(rq, sx);
    scm_gc_unprotect_object(sx);
  }
  void mcmclib_mh_q_guile_free(void* in_p) {
    /*scm_gc_unprotect_object((SCM) in_p);*/ //FIXME: find the correct way to
					     //release 'in_p'
  }
%}

%newobject mcmclib_mh_q_guile_alloc;
%inline %{
  mcmclib_mh_q* mcmclib_mh_q_guile_alloc(gsl_rng* r, SCM rq_dq) {
    scm_gc_protect_object(rq_dq);
    return mcmclib_mh_q_alloc(r,
			      mcmclib_mh_q_guile_rq,
			      mcmclib_mh_q_guile_dq,
			      (void*) rq_dq, mcmclib_mh_q_guile_free);
  }
%}

/**\brief Metropolis-Hastings proposal kernel*/
typedef struct mcmclib_mh_q_t {
  gsl_rng* r; /**< RNG used by the sampler function*/
  sampler_fun_t rq; /**< proposal sampler fun*/
  proposal_distr_fun_t dq; /**< proposal density fun*/
  void* gamma; /**< extra proposal kernel data*/
  free_fun_t free_gamma_fun; /**< optional gamma de-allocator fun*/
} mcmclib_mh_q;

%extend mcmclib_mh_q {
  ~mcmclib_mh_q() {
    mcmclib_mh_q_free($self);
  }
}

void mcmclib_mh_q_sample(mcmclib_mh_q* p, gsl_vector* x);
double mcmclib_mh_q_logd(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y);
double mcmclib_mh_q_ratio_offset(mcmclib_mh_q* p, gsl_vector* x, gsl_vector* y);
