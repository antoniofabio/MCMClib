%{
#include <amh.h>
#include <gauss_am.h>
#include <rapt.h>
#include <raptor.h>
%}

typedef void (*free_fun_t) (void* data);

typedef void (*mcmclib_amh_update_gamma_p) (void* p);

/**\brief Generic Adaptive Metropolis-Hastings sampler */
typedef struct {
  mcmclib_mh* mh;
  void* suff; /**< sufficient data accumulated up to current iteration*/
  mcmclib_amh_update_gamma_p update_gamma;
  int n; /**< current iteration number*/
} mcmclib_amh;

%extend mcmclib_amh {
  mcmclib_amh(mcmclib_mh* mh, void* suff,
	      mcmclib_amh_update_gamma_p update_gamma,
	      free_fun_t free_suff) {
    return mcmclib_amh_alloc(mh, suff, update_gamma, free_suff);
  }
  ~mcmclib_amh() {
    mcmclib_amh_free($self);
  }
}
int mcmclib_amh_update(mcmclib_amh* p);
void mcmclib_amh_reset(mcmclib_amh* p);

/*Concrete AMH samplers implementations*/

/*Gaussian AM*/
%newobject mcmclib_gauss_am_alloc;
mcmclib_amh* mcmclib_gauss_am_alloc(gsl_rng* r,
				    distrfun_p f, void* data,
				    gsl_vector* start_x,
				    const gsl_matrix* sigma_zero, int t0);
void mcmclib_gauss_am_set_sf(mcmclib_amh* p, double sf);

/*RAPT*/
typedef int (*region_fun_t) (void*, gsl_vector*);
%{
  int mcmclib_guile_region_fun(void* p, gsl_vector* x) {
    SCM sx = SWIG_NewPointerObj(x, SWIGTYPE_p_gsl_vector, 0);
    SCM ans = scm_call_1((SCM) p, sx);
    return scm_to_int(ans);
  }
%}

%typemap(in) (void* which_region_data, region_fun_t which_region) {
  scm_permanent_object($input); /*FIXME. Maybe solve it by proper use of '$owner'*/
  $1 = mcmclib_guile_region_fun;
  $2 = (void*) $input;
}

%newobject mcmclib_rapt_alloc;
mcmclib_amh* mcmclib_rapt_alloc(gsl_rng* r,
				distrfun_p f, void* data,
				gsl_vector* x,
				int t0,
				const gsl_matrix* sigma_whole,
				int K,
				gsl_matrix** sigma_local,
				region_fun_t which_region,
				void* which_region_data);

/*RAPTOR*/
%newobject mcmclib_raptor_alloc;
mcmclib_amh* mcmclib_raptor_alloc(gsl_rng* r,
				  distrfun_p f, void* data,
				  gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
				  gsl_vector* beta_hat,
				  gsl_vector** mu_hat,
				  gsl_matrix** Sigma_hat);
typedef double (*mcmclib_raptor_alpha_fun_t) (void* data, mcmclib_raptor_gamma*);
/** \brief RAPTOR sampler gamma values */
typedef struct {
  gsl_vector* beta_hat; /**< current mixture weights estimates*/
  gsl_vector** mu_hat; /**< current mixture means estimates*/
  gsl_matrix** Sigma_hat; /**< current mixture variances estimates*/

  mcmclib_mvnorm_lpdf** pik_hat; /**< single mixture components densities*/
  mcmclib_mixnorm_lpdf* pi_hat; /**< mixture density*/
} mcmclib_raptor_gamma;
%extend mcmclib_raptor_gamma {
  mcmclib_raptor_gamma(gsl_vector* beta_hat,
		       gsl_vector** mu_hat,
		       gsl_matrix** Sigma_hat) {
    return mcmclib_raptor_gamma_alloc(beta_hat, mu_hat, Sigma_hat);
  }
  ~mcmclib_raptor_gamma() {
    mcmclib_raptor_gamma_free($self);
  }
}
void mcmclib_raptor_set_sf(mcmclib_amh* p, double sf);
void mcmclib_raptor_set_sf_global(mcmclib_amh* p, double sf);
void mcmclib_raptor_set_sf_local(mcmclib_amh* p, double sf);
void mcmclib_raptor_set_alpha(mcmclib_amh* p, double alpha);

%{
#define RAPT_GAMMA(p) ((mcmclib_rapt_gamma*) (p)->mh->q->gamma)
#define RAPTOR_GAMMA(p) ((mcmclib_raptor_gamma*) RAPT_GAMMA(p)->which_region_data)
#define RAPTOR_SUFF(p) ((mcmclib_raptor_suff*) (p)->suff)
%}
%inline{
  mcmclib_raptor_gamma* mcmclib_raptor_gamma_get(mcmclib_amh* p) {
    return RAPTOR_GAMMA(p);
  }
}

%include "at7.i"
