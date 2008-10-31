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

  return ans;
}

void mcmclib_rapt_free(mcmclib_rapt* p) {
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

/*parametrized log-density of the (mixture) proposal function.
  To be used for computing M-H ratio correctly when doing the metropolis step
 */
static double q(void* data, gsl_vector* x, gsl_vector* y) {
  mcmclib_rapt* p = (mcmclib_rapt*) data;
  int region_x = p->which_region(x, p->which_region_data);
  mcmclib_mvnorm_lpdf* distr_obj =
    mcmclib_mvnorm_lpdf_alloc(x, (p->sigma_local[region_x])->data);
  double ans = mcmclib_mvnorm_lpdf_compute(distr_obj, y);
  mcmclib_mvnorm_lpdf_free(distr_obj);
  return ans;
}

int mcmclib_rapt_update(mcmclib_rapt* p) {
  gsl_rng* r = p->r;
  int *t = &(p->t);
  int t0 = p->t0;
  gsl_vector* old = p->old;
  gsl_vector* x = p->current_x;
  distrfun_p logdistr = p->logdistr;
  void* logdistr_data = p->logdistr_data;
  gsl_matrix* lambda = p->lambda;
  int K = p->K;
  gsl_matrix** sigma_local = p->sigma_local;
  gsl_matrix* sigma_whole = p->sigma_whole;
  gsl_matrix** variances = p->variances;
  gsl_vector** means = p->means;
  gsl_vector* n = p->n;
  gsl_matrix* visits = p->visits;
  region_fun_t which_region = p->which_region;
  void* which_region_data = p->which_region_data;
  gsl_matrix* jd = p->jd;

  gsl_vector_memcpy(old, x); /*save old state*/
  int which_region_old = which_region(old, which_region_data);
  gsl_vector_view lambda_vw = gsl_matrix_row(lambda, which_region_old);
  gsl_vector* lambda_p = &(lambda_vw.vector);

  /*update current chain value*/
  int which_proposal = sample(r, lambda_p); /*sample an integer between 0 and K, with given probabilities*/
  mcmclib_mvnorm(r,
		 (which_proposal < K) ? sigma_local[which_proposal] : sigma_whole,
		 x);
  gsl_vector_add(x, old);
  int accepted = mcmclib_metropolis_generic_step(r, old, x, logdistr, logdistr_data, q, p);
  int which_region_x = accepted ? which_region(x, which_region_data) : which_region_old;
  int k = which_region_x;

  /*update means and variances*/
  int fake_n = gsl_vector_get(n, k);
  mcmclib_covariance_update(variances[k], means[k], &fake_n, x);
  fake_n = (*t);
  mcmclib_covariance_update(p->global_variance, p->global_mean, &fake_n, x);

  /*update visits counts*/
  (*t)++;
  gsl_vector_set(n, which_region_x, gsl_vector_get(n, which_region_x) + 1);
  gsl_matrix_set(visits, which_region_x, which_proposal,
		 gsl_matrix_get(visits, which_region_x, which_proposal) + 1);

  /*if newreg == oldreg, update jumping distances*/
  if(which_region_old == which_region_x) {
    gsl_vector_sub(old, x);
    double newjd = 0.0;
    for(int i=0; i< x->size; i++)
      newjd += (gsl_vector_get(old, i) * gsl_vector_get(old, i));
    double newvisits = gsl_matrix_get(visits, k, which_proposal);
    gsl_matrix_set(jd, k, which_proposal,
		   (gsl_matrix_get(jd, k, which_proposal) *
		    (newvisits - 1) + newjd) / newvisits);
  }

  /*adaptation code*/
  if((*t) > t0) {
    /*if jd changed, update lambda weights*/
    int all_visited = 1;
    for(int j=0; j<=K; j++)
      if(gsl_matrix_get(visits, k, j)==0) {
	all_visited = 0;
	break;
      }
    if((which_region_old == which_region_x) && all_visited) {
      double sumjd = 0.0;
      for(int j=0; j<=K; j++)
	sumjd += gsl_matrix_get(jd, k, j);
      for(int j=0; j<=K; j++)
	gsl_matrix_set(lambda, k, j,
		       gsl_matrix_get(jd, k, j) / sumjd);
    }

    /*update local and global proposals covariance matrices*/
    gsl_matrix_memcpy(sigma_local[which_region_x], variances[which_region_x]);
    gsl_matrix_add(sigma_local[which_region_x], p->Sigma_eps);
    gsl_matrix_scale(sigma_local[which_region_x], 2.38 * 2.38 / ((double) x->size));
    gsl_matrix_memcpy(sigma_whole, p->global_variance);
    gsl_matrix_add(sigma_whole, p->Sigma_eps);
    gsl_matrix_scale(sigma_whole, 2.38 * 2.38 / ((double) x->size));
  }

  return 1;
}
