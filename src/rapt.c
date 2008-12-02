#include <gsl/gsl_math.h>
#include "metropolis.h"
#include "rapt.h"
#include "mvnorm.h"
#include "vector_stats.h"

static void rapt_init(mcmclib_rapt*);

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
  ans->means = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  ans->variances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    ans->means[k] = gsl_vector_alloc(dim);
    ans->variances[k] = gsl_matrix_alloc(dim, dim);
  }
  ans->global_mean = gsl_vector_alloc(dim);
  ans->global_variance = gsl_matrix_alloc(dim, dim);
  ans->n = gsl_vector_alloc(K);
  ans->lambda = gsl_matrix_alloc(K, K+1);
  ans->Sigma_eps = gsl_matrix_alloc(dim, dim);

  /*alloc extra data for mixture proposal density computation*/
  ans->q_mean = gsl_vector_alloc(dim);
  ans->q_k = (mcmclib_mvnorm_lpdf**) malloc((K+1) * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++)
      ans->q_k[k] = mcmclib_mvnorm_lpdf_alloc(ans->q_mean, ans->sigma_local[k]->data);
  ans->q_k[K] = mcmclib_mvnorm_lpdf_alloc(ans->q_mean, ans->sigma_whole->data);

  ans->workspace = gsl_vector_alloc(dim);

  rapt_init(ans);

  return ans;
}

static void rapt_init(mcmclib_rapt* p) {
  p->accepted = 1;
  p->t = 0;
  gsl_vector_set_all(p->global_mean, 0.0);
  gsl_matrix_set_all(p->global_variance, 0.0);
  gsl_vector_set_all(p->n, 0.0);
  gsl_matrix_set_all(p->lambda, 0.0);
  for(int region=0; region < p->K; region++) {
    gsl_matrix_set(p->lambda, region, region, 0.5);
    gsl_matrix_set(p->lambda, region, p->K, 0.5);
  }
  gsl_matrix_set_identity(p->Sigma_eps);
  gsl_matrix_scale(p->Sigma_eps, 0.001);

  p->which_region_x = p->which_region(p->current_x, p->which_region_data);
}

void mcmclib_rapt_free(mcmclib_rapt* p) {
  /*internal data free*/
  gsl_vector_free(p->workspace);
  gsl_matrix_free(p->Sigma_eps);
  gsl_matrix_free(p->lambda);
  gsl_vector_free(p->n);
  for(int k=0; k< p->K; k++) {
    gsl_vector_free(p->means[k]);
    gsl_matrix_free(p->variances[k]);
  }
  gsl_vector_free(p->global_mean);
  gsl_matrix_free(p->global_variance);
  free(p->means);
  free(p->variances);

  /*free proposal mixture density extra data*/
  gsl_vector_free(p->q_mean);
  for(int k=0; k <= p->K; k++)
    mcmclib_mvnorm_lpdf_free(p->q_k[k]);
  free(p->q_k);

  gsl_matrix_free(p->sigma_whole);
  for(int k=0; k < p->K; k++)
    gsl_matrix_free(p->sigma_local[k]);
  free(p->sigma_local);

  gsl_vector_free(p->old);
  free(p);
}

/**
Follows a bunch of prototypes of internal funcs. used by rapt_update
*/
/*sample a discrete value from the discrete distribution with probs. 'probs'*/
static int sample(gsl_rng* r, gsl_vector* probs);
/*log-density of the (mixture) proposal kernel*/
static double rapt_q(void* data, gsl_vector* x, gsl_vector* y);
/*does the metropolis step*/
static void rapt_update_current_value(mcmclib_rapt* p);
/*update current means and variances values basing on last chain step*/
static void rapt_update_means_variances(mcmclib_rapt* p);

int mcmclib_rapt_update(mcmclib_rapt* p) {
  /*update current chain value*/
  rapt_update_current_value(p);
  /*update means and variances*/
  rapt_update_means_variances(p);

  return 1;
}

static void rapt_update_current_value(mcmclib_rapt* p) {
  /*save old state, old region*/
  gsl_vector_memcpy(p->old, p->current_x);
  p->which_region_old = p->which_region_x;

  gsl_vector_view lambda_vw = gsl_matrix_row(p->lambda, p->which_region_old);
  gsl_vector* lambda_p = &(lambda_vw.vector);

  /*sample an integer between 0 and K, with given probabilities:*/
  p->which_proposal = sample(p->r, lambda_p);

  /*sample an error term from the corresponding prob. distrib.*/
  mcmclib_mvnorm(p->r, (p->which_proposal < p->K) ?
		 p->sigma_local[p->which_proposal] : p->sigma_whole,
		 p->current_x);
  gsl_vector_add(p->current_x, p->old);
  p->accepted = mcmclib_metropolis_generic_step(p->r,
						p->old,
						p->current_x,
						p->logdistr, p->logdistr_data,
						rapt_q, p);
  p->which_region_x = p->accepted ? p->which_region(p->current_x, p->which_region_data) : p->which_region_old;

  (p->t)++;
}

/*parametrized log-density of the (mixture) proposal function.
  To be used for computing M-H ratio correctly when doing the metropolis step
 */
static double rapt_q(void* data, gsl_vector* x, gsl_vector* y) {
  mcmclib_rapt* p = (mcmclib_rapt*) data;
  int region_x = p->which_region(x, p->which_region_data);
  gsl_vector_view lambda_vw = gsl_matrix_row(p->lambda, region_x);
  gsl_vector* lambda_p = &(lambda_vw.vector);

  double ans = 0.0;
  gsl_vector_memcpy(p->q_mean, x);
  for(int k=0; k <= p->K; k++)
    ans += exp(mcmclib_mvnorm_lpdf_compute(p->q_k[k], y)) * gsl_vector_get(lambda_p, k);

  return log(ans);
}

static int sample(gsl_rng* r, gsl_vector* probs) {
  int K = probs->size;
  double cum_sum = 0.0;
  double who = gsl_rng_uniform(r);
  for(int which=0; which<K; which++) {
    if(who < (cum_sum += gsl_vector_get(probs, which)))
      return(which);
  }
  return(K-1);
}

/*dest = alpha * (A + B) */
static void matrix_addscale(gsl_matrix* dest,
			     gsl_matrix* A, gsl_matrix* B, double alpha) {
  gsl_matrix_memcpy(dest, A);
  gsl_matrix_add(dest, B);
  gsl_matrix_scale(dest, alpha);
}

void mcmclib_rapt_update_proposals_custom(mcmclib_rapt* p,
					  gsl_matrix** variances,
					  gsl_matrix* global_variance) {
  double sf = 2.38 * 2.38 / ((double) p->old->size);
  for(int k=0; k< p->K; k++)
    matrix_addscale(p->sigma_local[k],
			  variances[k], p->Sigma_eps, sf);
  matrix_addscale(p->sigma_whole,
			global_variance, p->Sigma_eps, sf);
}

void mcmclib_rapt_update_proposals(mcmclib_rapt* p) {
  if((p->t) <= p->t0)
    return;
  mcmclib_rapt_update_proposals_custom(p, p->variances, p->global_variance);
}

static void rapt_update_means_variances(mcmclib_rapt* p) {
  int k = p->which_region_x;
  int fake_n = gsl_vector_get(p->n, k);
  mcmclib_covariance_update(p->variances[k], p->means[k], &fake_n, p->current_x);
  fake_n = p->t;
  mcmclib_covariance_update(p->global_variance, p->global_mean, &fake_n, p->current_x);
}
