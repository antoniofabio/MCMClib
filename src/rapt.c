#include <gsl/gsl_math.h>
#include "metropolis.h"
#include "rapt.h"
#include "mvnorm.h"

mcmclib_rapt* mcmclib_rapt_alloc(
				 gsl_rng* r,
				 distrfun_p logdistr, void* logdistr_data,
				 gsl_vector* x,
				 int t0,
				 const gsl_matrix* sigma_whole,
				 int K,
				 gsl_matrix** sigma_local,
				 region_fun_t which_region,
				 void* which_region_data){
  int dim = x->size;
  mcmclib_rapt* ans = (mcmclib_rapt*) malloc(sizeof(mcmclib_rapt));
  ans->r = r;
  ans->logdistr = logdistr;
  ans->logdistr_data = logdistr_data;
  ans->current_x = x;
  ans->old = gsl_vector_alloc(dim);

  ans->t0 = t0;
  ans->sigma_whole = gsl_matrix_alloc(dim, dim);
  gsl_matrix_memcpy(ans->sigma_whole, sigma_whole);
  ans->K = K;
  ans->sigma_local = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    ans->sigma_local[k] = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(ans->sigma_local[k], sigma_local[k]);
  }
  ans->which_region = which_region;
  ans->which_region_data = which_region_data;

  /*internal data alloc*/
  ans->t = 0;
  ans->whole_variance = gsl_matrix_alloc(dim, dim);
  ans->means = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  ans->variances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    ans->means[k] = gsl_vector_alloc(dim);
    ans->variances[k] = gsl_matrix_alloc(dim, dim);
  }
  ans->visits = gsl_matrix_alloc(K, K+1);
  gsl_matrix_set_all(ans->visits, 0.0);
  ans->jd = gsl_vector_alloc(K);
  gsl_vector_set_all(ans->jd, 0.0);
  ans->n = gsl_vector_alloc(K);
  gsl_vector_set_all(ans->n, 0.0);
  ans->lambda = gsl_vector_alloc(K+1);
  gsl_vector_set_all(ans->lambda, 1.0 / (double) (K+1.0));

  return ans;
}

void mcmclib_rapt_free(mcmclib_rapt* p) {
  /*internal data free*/
  gsl_vector_free(p->lambda);
  gsl_matrix_free(p->visits);
  gsl_vector_free(p->jd);
  gsl_vector_free(p->n);
  gsl_matrix_free(p->whole_variance);
  for(int k=0; k< p->K; k++) {
    gsl_vector_free(p->means[k]);
    gsl_matrix_free(p->variances[k]);
  }
  free(p->means);
  free(p->variances);

  gsl_matrix_free(p->sigma_whole);
  for(int k=0; k < p->K; k++)
    gsl_matrix_free(p->sigma_local[k]);
  free(p->sigma_local);

  gsl_vector_free(p->old);
  free(p);
}

/*sample a discrete value from the discrete distribution with probs. 'probs'*/
static int sample(gsl_rng* r, gsl_vector* probs) {
  int K = probs->size;
  double cum_sum = 0.0;
  double who = gsl_rng_uniform(r);
  for(int which=0; which<K; which++) {
    cum_sum += gsl_vector_get(probs, which);
    if(who > cum_sum)
      return(which);
  }
  return(K-1);
}

/*TODO*/
int mcmclib_rapt_update(mcmclib_rapt* p) {
  gsl_rng* r = p->r;
  int *t = &(p->t);
  int t0 = p->t0;
  gsl_vector* x = p->current_x;
  distrfun_p logdistr = p->logdistr;
  void* logdistr_data = p->logdistr_data;
  gsl_vector* lambda = p->lambda;
  int K = p->K;
  gsl_matrix** sigma_local = p->sigma_local;
  gsl_matrix* sigma_whole = p->sigma_whole;

  /*step 1: update current state*/
  gsl_vector_memcpy(p->old, x); /*save old state*/
  int which_proposal = sample(r, lambda); /*sample an integer between 0 and K, with given probabilities*/
  mcmclib_mvnorm(r,
		 (which_proposal < K) ? sigma_local[which_proposal] : sigma_whole,
		 x);
  /*TODO: compute (correctly...) Metropolis ratio*/

  /*step 2: update means and variances*/
  /*step 3: update proposal covariance matrices*/
  /*step 4: update weights*/

  return 1;
}
