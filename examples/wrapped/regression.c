/**Wrapped gaussian regression example*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mh.h>
#include <gauss_am.h>

#define S 100 /*number of points*/
#define OUTPUT_FILE "chain.csv"
#define N 150000 /*chain length*/
#define THIN 10 /*thinning interval*/
#define K (3 + S) /*total number of parameters*/
#define V0 1.0 /*step size*/

#define BETA0_MU 0.0
#define BETA0_SIGMA 1.0
#define BETA1_MU 0.0
#define BETA1_SIGMA 1.0
#define SIGMA2_A 102.0
#define SIGMA2_B 1010.0

gsl_vector* y; /*observed data*/

gsl_vector* X; /*regressor*/

gsl_vector* x_theta; /*chain: intercept, slope, variance*/
gsl_vector* x_k;     /*chain: wrapping coefficients*/

void sampk(void* data, gsl_vector* out) {
  mcmclib_mh_q* q = data;
  for(int s=0; s<S; s++) {
    int coin = gsl_rng_uniform(q->r) <= 0.1;
    if(coin) {
      double ks = round(gsl_ran_flat(q->r, -1.5, 1.5));
      gsl_vector_set(out, s, gsl_vector_get(out, s) + ks);
    }
  }
}
double sampk_d(void* data, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

/*log-likelihood function*/
double loglik(gsl_vector* theta, gsl_vector* kv) {
  double beta0 = gsl_vector_get(theta, 0);
  double beta1 = gsl_vector_get(theta, 1);
  double sigma = sqrt(exp(gsl_vector_get(theta, 2)));
  double* k = kv->data;
  double ans = 0.0;
  if(sigma <= 0)
    return log(0.0);
  for(int i=0; i<S; i++) {
    double xi = gsl_vector_get(y, i);
    double mi = beta0 + beta1 * gsl_vector_get(X, i);
    ans += log(gsl_ran_gaussian_pdf(xi - M_PI + 2.0 * k[i] * M_PI - mi, sigma));
  }
  return ans;
}

double linvGamma(double x, double a, double b) {
  return a * log(b) - (a + 1.0) * log(x) - (b / x) - gsl_sf_lngamma(a);
}

/*target (beta, sigma) distribution*/
double post_theta(void* ignore, gsl_vector* x) {
  return loglik(x, x_k)
    + log(gsl_ran_gaussian_pdf(gsl_vector_get(x, 0) - BETA0_MU, sqrt(BETA0_SIGMA)))
    + log(gsl_ran_gaussian_pdf(gsl_vector_get(x, 1) - BETA1_MU, sqrt(BETA1_SIGMA)))
    + linvGamma(exp(gsl_vector_get(x, 2)), SIGMA2_A, SIGMA2_B);
}

/*target k distribution*/
double post_k(void* ignore, gsl_vector* x) {
  return loglik(x_theta, x);
}

int main(int argc, char** argv) {
  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  /*set regressor data*/
  X = gsl_vector_alloc(S);
  for(int i=0; i<S; i++)
    gsl_vector_set(X, i, (double) i * 10.0 / (double) S);
  /*set observed data*/
  y = gsl_vector_alloc(S);
  for(int i=0; i<S; i++)
    gsl_vector_set(y, i, gsl_ran_flat(r, -M_PI, M_PI));

  /*set starting chain values*/
  x_theta = gsl_vector_alloc(3);
  gsl_vector_set(x_theta, 0, 0.0);
  gsl_vector_set(x_theta, 1, 0.0);
  gsl_vector_set(x_theta, 2, log(1000));
  x_k = gsl_vector_alloc(S);
  gsl_vector_set_all(x_k, 0.0);

  /*alloc a new sampler for theta*/
  int theta_d = x_theta->size;
  gsl_matrix* s_theta_Sigma0 = gsl_matrix_alloc(theta_d, theta_d);
  gsl_matrix_set_identity(s_theta_Sigma0);
  gsl_matrix_scale(s_theta_Sigma0, V0);
  mcmclib_amh* s_theta = mcmclib_gauss_am_alloc(r, post_theta, NULL, x_theta, s_theta_Sigma0, 10000);

  /*alloc a new sampler for k*/
  mcmclib_mh_q* q_k = mcmclib_mh_q_alloc(r, sampk, NULL, sampk_d, NULL, NULL);
  mcmclib_mh* s_k = mcmclib_mh_alloc(r, post_k, NULL, q_k, x_k);

  /*open output file*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  /*print out csv header*/
  fprintf(out, "beta0, beta1, sigma");
  for(int s=0; s<S; s++)
    fprintf(out, ", k%d", s+1);
  fprintf(out, "\n");

  /*main MCMC loop*/
  for(int i=0; i<N; i++) {
    mcmclib_amh_update(s_theta);
    mcmclib_mh_update(s_k);
    if(((i+1) % THIN) == 0) {
      for(int j=0; j<theta_d; j++)
	fprintf(out, "%f, ", gsl_vector_get(x_theta, j));
      for(int j=0; j<(S-1); j++)
	fprintf(out, "%f, ", gsl_vector_get(x_k, j));
      fprintf(out, "%f\n", gsl_vector_get(x_k, S-1));
    }
  }

  fclose(out);
  gsl_rng_free(r);
  gsl_vector_free(x_theta);
  gsl_vector_free(x_k);
  gsl_vector_free(y);
  gsl_vector_free(X);
  mcmclib_gauss_am_free(s_theta);
  mcmclib_mh_free(s_k);
  mcmclib_mh_q_free(q_k);
}
