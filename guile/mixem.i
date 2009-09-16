%{
#include <mixem_online.h>
%}

/**On-Line EM for gaussian mixtures*/
typedef struct {
  mcmclib_mixolem_suff* gamma; /**< pointer to current parameters estimates*/
} mcmclib_mixem_online;

%extend mcmclib_mixem_online {
  mcmclib_mixem_online(gsl_vector** mu,
		       gsl_matrix** Sigma,
		       gsl_vector* beta,
		       double eta_eps,
		       int n0) {
    return mcmclib_mixem_online_alloc(mu, Sigma, beta, eta_eps, n0);
  }
  ~mcmclib_mixem_online() {
    mcmclib_mixem_online_free($self);
  }
}

void mcmclib_mixem_online_update(mcmclib_mixem_online* p, gsl_vector* y);
