#include <gsl/gsl_math.h>
#include "rapt_q.h"
#include "mvnorm.h"

void mcmclib_rapt_gamma_set_alpha(mcmclib_rapt_gamma* p, double alpha) {
  gsl_matrix_set_all(p->lambda, 0.0);
  for(int region=0; region < p->K; region++) {
    gsl_matrix_set(p->lambda, region, region, 1.0 - alpha);
    gsl_matrix_set(p->lambda, region, p->K, alpha);
  }
}

mcmclib_rapt_gamma* mcmclib_rapt_gamma_alloc(const gsl_matrix* sigma_whole,
					     int K,
					     gsl_matrix** sigma_local,
					     region_fun_t which_region,
					     void* which_region_data) {
  int dim = sigma_whole->size1;
  mcmclib_rapt_gamma* a = (mcmclib_rapt_gamma*) malloc(sizeof(mcmclib_rapt_gamma));
  a->sigma_whole = gsl_matrix_alloc(dim, dim);
  gsl_matrix_memcpy(a->sigma_whole, sigma_whole);
  a->K = K;
  a->sigma_local = (gsl_matrix**) malloc(K * sizeof(gsl_matrix*));
  for(int k=0; k<K; k++) {
    a->sigma_local[k] = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(a->sigma_local[k], sigma_local[k]);
  }
  a->which_region = which_region;
  a->which_region_data = which_region_data;

  a->lambda = gsl_matrix_alloc(K, K+1);
  mcmclib_rapt_gamma_set_alpha(a, 0.5);

  a->q_mean = gsl_vector_alloc(dim);
  a->q_k = (mcmclib_mvnorm_lpdf**) malloc((K+1) * sizeof(mcmclib_mvnorm_lpdf*));
  for(int k=0; k<K; k++)
      a->q_k[k] = mcmclib_mvnorm_lpdf_alloc(a->q_mean, a->sigma_local[k]->data);
  a->q_k[K] = mcmclib_mvnorm_lpdf_alloc(a->q_mean, a->sigma_whole->data);
  a->workspace = gsl_vector_alloc(dim);
  return a;
}

void mcmclib_rapt_gamma_free(mcmclib_rapt_gamma* p) {
  gsl_vector_free(p->workspace);
  gsl_matrix_free(p->lambda);

  gsl_vector_free(p->q_mean);
  for(int k=0; k <= p->K; k++)
    mcmclib_mvnorm_lpdf_free(p->q_k[k]);
  free(p->q_k);

  gsl_matrix_free(p->sigma_whole);
  for(int k=0; k < p->K; k++)
    gsl_matrix_free(p->sigma_local[k]);
  free(p->sigma_local);

  free(p);
}

void mcmclib_rapt_q_sample(void* in_p, gsl_vector* x);
double mcmclib_rapt_q_d(void* gamma, gsl_vector* x, gsl_vector* y);

mcmclib_mh_q* mcmclib_rapt_q_alloc(gsl_rng* r,
				   distrfun_p logdistr, void* logdistr_data,
				   const gsl_matrix* sigma_whole,
				   int K,
				   gsl_matrix** sigma_local,
				   region_fun_t which_region,
				   void* which_region_data) {
  mcmclib_rapt_gamma* gamma = mcmclib_rapt_gamma_alloc(sigma_whole,
						       K, sigma_local,
						       which_region, which_region_data);
  return mcmclib_mh_q_alloc(r, mcmclib_rapt_q_sample, gamma,
			    mcmclib_rapt_q_d, gamma,
			    gamma);
}

void mcmclib_rapt_q_free(mcmclib_mh_q* p) {
  mcmclib_rapt_gamma_free(p->gamma);
  mcmclib_mh_q_free(p);
}

/*sample a discrete value from the discrete distribution with probs. 'probs'*/
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

void mcmclib_rapt_q_sample(void* in_p, gsl_vector* x) {
  mcmclib_mh_q* q = (mcmclib_mh_q*) in_p;
  mcmclib_rapt_gamma* g = (mcmclib_rapt_gamma*) q->gamma;

  g->which_region_old = g->which_region_x = g->which_region(x, g->which_region_data);
  gsl_vector_memcpy(g->workspace, x);

  gsl_vector_view lambda_vw = gsl_matrix_row(g->lambda, g->which_region_old);
  gsl_vector* lambda_p = &(lambda_vw.vector);

  /*sample an integer between 0 and K, with given probabilities:*/
  g->which_proposal = sample(q->r, lambda_p);

  /*sample an error term from the corresponding prob. distrib.*/
  mcmclib_mvnorm(q->r, (g->which_proposal < g->K) ?
		 g->sigma_local[g->which_proposal] : g->sigma_whole,
		 x);
  gsl_vector_add(x, g->workspace);
}

/*parametrized log-density of the (mixture) proposal function.
  To be used for computing M-H ratio correctly when doing the metropolis step
 */
double mcmclib_rapt_q_d(void* gamma, gsl_vector* x, gsl_vector* y) {
  mcmclib_rapt_gamma* p = (mcmclib_rapt_gamma*) gamma;
  int region_x = p->which_region(x, p->which_region_data);
  gsl_vector_view lambda_vw = gsl_matrix_row(p->lambda, region_x);
  gsl_vector* lambda_p = &(lambda_vw.vector);

  double ans = 0.0;
  gsl_vector_memcpy(p->q_mean, x);
  for(int k=0; k <= p->K; k++)
    ans += exp(mcmclib_mvnorm_lpdf_compute(p->q_k[k], y)) * gsl_vector_get(lambda_p, k);

  return log(ans);
}
