%{
#include <amh.h>
#include <gauss_am.h>
#include <raptor.h>
%}

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
	      mcmclib_amh_update_gamma_p update_gamma) {
    return mcmclib_amh_alloc(mh, suff, update_gamma);
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
%delobject mcmclib_gauss_am_free;
mcmclib_amh* mcmclib_gauss_am_alloc(gsl_rng* r,
				    distrfun_p logdistr, void* logdistr_data,
				    gsl_vector* start_x,
				    const gsl_matrix* sigma_zero, int t0);
void mcmclib_gauss_am_free(mcmclib_amh* p);
void mcmclib_gauss_am_set_sf(mcmclib_amh* p, double sf);

/*RAPTOR*/
mcmclib_amh* mcmclib_raptor_alloc(gsl_rng* r,
				  distrfun_p logdistr, void* logdistr_data,
				  gsl_vector* x, int t0, gsl_matrix* Sigma_zero,
				  gsl_vector* beta_hat,
				  gsl_vector** mu_hat,
				  gsl_matrix** Sigma_hat);
void mcmclib_raptor_free(mcmclib_amh* p);
