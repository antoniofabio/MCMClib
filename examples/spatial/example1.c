/**Spatial model example*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <spatial.h>
#include <gauss_rw.h>

#define S 27 /*number of points*/
#define INPUT_FILE "data.dat"
#define OUTPUT_FILE "chain.csv"
#define N 150000 /*chain length*/
#define THIN 10 /*thinning interval*/
#define V0 0.1 /*step size*/
#define K (3 + S) /*total number of parameters*/

#define MU_MU0 0.5
#define MU_SIGMA2 0.1

gsl_vector* theta; /*spatial covariance parameters vector*/
gsl_vector_view rho_v; /*range*/
gsl_vector_view sigma_v; /*sill*/
gsl_vector_view tausq_v; /*nugget*/
gsl_vector_view mu_v; /*mean*/

gsl_vector* y_obs; /*observed data*/

mcmclib_spatial_lpdf* lik; /*data likelihood object*/

/*log-likelihood function*/
double loglik(gsl_vector* in_theta) {
  gsl_vector_memcpy(theta, in_theta);
  return mcmclib_spatial_lpdf_compute(lik, y_obs);
}

double linvGamma(double x, double a, double b) {
  return a * log(b) - (a + 1.0) * log(x) - (b / x) - gsl_sf_lngamma(a);
}

double lmu_prior(gsl_vector* x, double mu0, double sigma0) {
  int n = x->size;
  double ans = 0.0;
  for(int s=0; s<n; s++)
    ans += log(gsl_ran_gaussian_pdf(gsl_vector_get(x, s) - mu0, sqrt(sigma0)));
  return ans;
}

/*target distribution*/
double target_logdensity(void* ignore, gsl_vector* x) {
  double r = x->data[0];
  double s = x->data[1];
  double t = x->data[2];
  gsl_vector_view mu_v = gsl_vector_subvector(x, 3, S);
  gsl_vector* mu = &(mu_v.vector);
  if((r<0) || (s<0) || (t<0) )
    return log(0.0);
  return loglik(x)
    + gsl_ran_gamma_pdf(r, 3.0, 40.0)
    + linvGamma(s, 8.0, 9.0)
    + gsl_ran_gamma_pdf(t, 0.1, 0.1)
    + lmu_prior(mu, MU_MU0, MU_SIGMA2);
}

int main(int argc, char** argv) {
  /*read input data*/
  gsl_matrix* data = gsl_matrix_alloc(S, 3);
  gsl_matrix_view xy_v = gsl_matrix_submatrix(data, 0, 0, S, 2);
  gsl_matrix* xy = &(xy_v.matrix);
  gsl_vector_view y_obs_v = gsl_matrix_column(data, 2);
  y_obs = &(y_obs_v.vector);
  gsl_matrix* D = gsl_matrix_alloc(S, S);
  mcmclib_spatial_distances(D, xy);

  /*set model parameters*/
  theta = gsl_vector_alloc(K);
  rho_v = gsl_vector_subvector(theta, 0, 1);
  sigma_v = gsl_vector_subvector(theta, 1, 1);
  tausq_v = gsl_vector_subvector(theta, 2, 1);
  mu_v = gsl_vector_subvector(theta, 3, S);
  lik = mcmclib_spatial_lpdf_alloc(&(mu_v.vector), &(rho_v.vector),
				   &(sigma_v.vector), &(tausq_v.vector), D);

  /*set starting chain values*/
  gsl_vector* x = gsl_vector_alloc(K);
  gsl_vector_set_all(x, 0.5);

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*alloc a new sampler*/
  mcmclib_mh* sampler = mcmclib_gauss_rw_alloc(r, target_logdensity, NULL, x, V0);

  /*open output file*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  /*print out csv header*/
  fprintf(out, "rho, sigma, tausq");
  for(int s=0; s<S; s++)
    fprintf(out, ", mu%d", s+1);
  fprintf(out, "\n");

  /*main MCMC loop*/
  for(int i=0; i<N; i++) {
    mcmclib_mh_update(sampler);
    if(((i+1) % THIN) == 0) {
      for(int j=0; j<(K-1); j++)
	fprintf(out, "%f, ", gsl_vector_get(x, j));
      fprintf(out, "%f\n", gsl_vector_get(x, K-1));
    }
  }

  fclose(out);
  gsl_rng_free(r);
  gsl_vector_free(x);
  mcmclib_gauss_rw_free(sampler);
  mcmclib_spatial_lpdf_free(lik);
  gsl_matrix_free(D);
  gsl_vector_free(theta);
  gsl_matrix_free(data);
}
