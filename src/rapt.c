#include <gsl/gsl_math.h>
#include "metropolis.h"
#include "rapt.h"
#include "mvnorm.h"
#include "vector_stats.h"

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

  ans->accepted = 1;
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
  ans->means = (gsl_vector**) malloc(K * sizeof(gsl_vector*));
  ans->variances = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    ans->means[k] = gsl_vector_alloc(dim);
    ans->variances[k] = gsl_matrix_alloc(dim, dim);
  }
  ans->global_mean = gsl_vector_alloc(dim);
  gsl_vector_set_all(ans->global_mean, 0.0);
  ans->global_variance = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_all(ans->global_variance, 0.0);
  ans->visits = gsl_matrix_alloc(K, K+1);
  gsl_matrix_set_all(ans->visits, 0.0);
  ans->jd = gsl_matrix_alloc(K, K+1);
  gsl_matrix_set_all(ans->jd, 0.0);
  ans->n = gsl_vector_alloc(K);
  gsl_vector_set_all(ans->n, 0.0);
  ans->lambda = gsl_matrix_alloc(K, K+1);
  gsl_matrix_set_all(ans->lambda, 1.0 / (double) (K+1.0));
  ans->Sigma_eps = gsl_matrix_alloc(dim, dim);
  gsl_matrix_set_identity(ans->Sigma_eps);
  gsl_matrix_scale(ans->Sigma_eps, 0.001);

  ans->ntries = gsl_vector_alloc(K+1);
  gsl_vector_set_all(ans->ntries, 0.0);
  ans->workspace = gsl_vector_alloc(dim);

  return ans;
}

void mcmclib_rapt_free(mcmclib_rapt* p) {
  /*extra data free*/
  gsl_vector_free(p->workspace);
  gsl_vector_free(p->ntries);

  /*internal data free*/
  gsl_matrix_free(p->Sigma_eps);
  gsl_matrix_free(p->lambda);
  gsl_matrix_free(p->visits);
  gsl_matrix_free(p->jd);
  gsl_vector_free(p->n);
  for(int k=0; k< p->K; k++) {
    gsl_vector_free(p->means[k]);
    gsl_matrix_free(p->variances[k]);
  }
  gsl_vector_free(p->global_mean);
  gsl_matrix_free(p->global_variance);
  free(p->means);
  free(p->variances);

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
/*update num. of trials from each proposal info relative to current point*/
static void rapt_update_ntries(mcmclib_rapt* p);
/*update visits counts info relative to each region, each proposal*/
static void rapt_update_visits_counts(mcmclib_rapt* p);
/*update mean jumping distance info relative to each region, each proposal*/
static void rapt_update_jumping_distances(mcmclib_rapt* p);
/***/

int mcmclib_rapt_update(mcmclib_rapt* p) {
  if(p->accepted == 1)
    gsl_vector_set_all(p->ntries, 0.0);

  /*update current chain value*/
  rapt_update_current_value(p);
  /*update means and variances*/
  rapt_update_means_variances(p);

  /*update ntries*/
  rapt_update_ntries(p);
  /*update visits counts*/
  rapt_update_visits_counts(p);
  /*update jumping distances*/
  rapt_update_jumping_distances(p);

  return 1;
}

void mcmclib_rapt_update_proposals(mcmclib_rapt* p) {
  if((p->t) > p->t0) {
    gsl_matrix_memcpy(p->sigma_local[p->which_region_x], p->variances[p->which_region_x]);
    gsl_matrix_add(p->sigma_local[p->which_region_x], p->Sigma_eps);
    gsl_matrix_scale(p->sigma_local[p->which_region_x], 2.38 * 2.38 / ((double) p->old->size));
    gsl_matrix_memcpy(p->sigma_whole, p->global_variance);
    gsl_matrix_add(p->sigma_whole, p->Sigma_eps);
    gsl_matrix_scale(p->sigma_whole, 2.38 * 2.38 / ((double) p->old->size));
  }
}

void mcmclib_rapt_update_lambda(mcmclib_rapt* p) {
  if(p->t <= p->t0)
    return;

  for(int k=0; k< p->K; k++) { /*for each region*/
    gsl_vector_view vrv = gsl_matrix_row(p->visits, k);
    gsl_vector* visits_k = &(vrv.vector);
    if(!gsl_vector_ispos(visits_k))
      continue;

    double beta = gsl_matrix_get(p->lambda, k, p->K);
    double sumjd = 0.0;
    for(int j=0; j< p->K; j++) /*get 'total' jumping distance*/
      sumjd += gsl_matrix_get(p->jd, k, j);
    for(int j=0; j< p->K; j++) /*for each proposal*/
      gsl_matrix_set(p->lambda, k, j,
		     (gsl_matrix_get(p->jd, k, j) / sumjd) * (1.0 - beta));
  }
}

static void rapt_update_ntries(mcmclib_rapt* p) {
  gsl_vector_set(p->ntries, p->which_proposal,
		 gsl_vector_get(p->ntries, p->which_proposal) + 1);
}

static void rapt_update_current_value(mcmclib_rapt* p) {
  /*save old state, old region*/
  gsl_vector_memcpy(p->old, p->current_x);
  p->which_region_old = p->which_region(p->old, p->which_region_data);

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
}

static void rapt_update_means_variances(mcmclib_rapt* p) {
  int k = p->which_region_x;
  int fake_n = gsl_vector_get(p->n, k);
  mcmclib_covariance_update(p->variances[k], p->means[k], &fake_n, p->current_x);
  fake_n = p->t;
  mcmclib_covariance_update(p->global_variance, p->global_mean, &fake_n, p->current_x);
}

static void rapt_update_visits_counts(mcmclib_rapt* p) {
  p->t++;
  gsl_vector_set(p->n, p->which_region_x,
		 gsl_vector_get(p->n, p->which_region_x) + 1);
  gsl_matrix_set(p->visits, p->which_region_x, p->which_proposal,
		 gsl_matrix_get(p->visits, p->which_region_x, p->which_proposal) + 1);

}

static void rapt_update_jumping_distances(mcmclib_rapt* p) {
  if(p->which_region_old != p->which_region_x)
    return;

  gsl_vector_memcpy(p->workspace, p->old);
  gsl_vector_sub(p->workspace, p->current_x);
  double newjd = 0.0;
  for(int i=0; i< p->current_x->size; i++)
    newjd += (gsl_vector_get(p->workspace, i) * gsl_vector_get(p->workspace, i));
  p->last_jd = newjd;
  int k = p->which_region_x;
  double newvisits = gsl_matrix_get(p->visits, k, p->which_proposal);
  gsl_matrix_set(p->jd, k, p->which_proposal,
		 (gsl_matrix_get(p->jd, k, p->which_proposal) *
		  (newvisits - 1) + newjd) / newvisits);
}

/*parametrized log-density of the (mixture) proposal function.
  To be used for computing M-H ratio correctly when doing the metropolis step
 */
static double rapt_q(void* data, gsl_vector* x, gsl_vector* y) {
  mcmclib_rapt* p = (mcmclib_rapt*) data;
  int region_x = p->which_region(x, p->which_region_data);
  mcmclib_mvnorm_lpdf* distr_obj =
    mcmclib_mvnorm_lpdf_alloc(x, (p->sigma_local[region_x])->data);
  double ans = mcmclib_mvnorm_lpdf_compute(distr_obj, y);
  mcmclib_mvnorm_lpdf_free(distr_obj);
  return ans;
}

static int sample(gsl_rng* r, gsl_vector* probs) {
  int K = probs->size;
  double cum_sum = 0.0;
  double who = gsl_rng_uniform(r);
  for(int which=0; which<K; which++) {
    cum_sum += gsl_vector_get(probs, which);
    if((who > cum_sum) &&
       (who <= (cum_sum + gsl_vector_get(probs, which+1))))
      return(which);
  }
  return(K-1);
}
