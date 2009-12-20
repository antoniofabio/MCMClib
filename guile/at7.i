%{
#include <at7.h>
%}

/*AT7 sampler*/
/** \brief AT7 gamma values */
typedef struct {
  gsl_vector* beta; /**< current mixture weights*/
  gsl_vector** mu; /**< current mixture means*/
  gsl_matrix** Sigma; /**< current mixture variances*/

  mcmclib_mixnorm_lpdf* pi; /**< fitted mixture density*/

  gsl_matrix** qVariances; /**< local proposal variances*/
  mcmclib_mvnorm_lpdf** qdk; /**< proposal density components*/

  gsl_matrix* Sigma_eps; /**< positive-definiteness correction factor*/
  gsl_vector* scaling_factors; /**< region-specific scaling factors*/
} mcmclib_at7_gamma;

/** alloc a new at7_gamma object. Input arguments are copied @internal */
mcmclib_at7_gamma* mcmclib_at7_gamma_alloc(const gsl_vector* beta_hat,
					   gsl_vector** mu_hat,
					   gsl_matrix** Sigma_hat);

/** \brief AT7 sufficient data */
typedef mcmclib_mixem_online* mcmclib_at7_suff;

%newobject mcmclib_at7_alloc;
mcmclib_amh* mcmclib_at7_alloc(gsl_rng* r,
			       distrfun_p f, void* data,
			       gsl_vector* x, int t0,
			       const gsl_vector* beta_hat,
			       gsl_vector** mu_hat,
			       gsl_matrix** Sigma_hat);
void mcmclib_at7_set_sf(mcmclib_amh* p, const gsl_vector* sf);
void mcmclib_at7_set_sf_all(mcmclib_amh* p, double all);

%inline{
  mcmclib_at7_gamma* mcmclib_at7_gamma_get(mcmclib_amh* p) {
    return ((mcmclib_at7_gamma*)((mcmclib_amh*) p)->mh->q->gamma);
  }
}
