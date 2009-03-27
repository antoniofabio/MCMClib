/**Spatial model example*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <spatial.h>
#include <gauss_rw.h>

#define S 9 /*number of points*/
#define OUTPUT_FILE "data.csv"
#define N 100000 /*chain length*/
#define V0 0.1 /*step size*/

gsl_vector* theta; /*spatial covariance parameters vector*/
gsl_vector_view rho_v; /*range*/
gsl_vector_view sigma_v; /*sill*/
gsl_vector_view tausq_v; /*nugget*/

gsl_vector* y_obs; /*observed data*/

mcmclib_spatial_lpdf* lik; /*data likelihood object*/

/*log-likelihood function*/
double loglik(gsl_vector* in_theta) {
  gsl_vector_memcpy(theta, in_theta);
  return mcmclib_spatial_lpdf_compute(lik, y_obs);
}

/*target distribution*/
double target_logdensity(void* ignore, gsl_vector* x) {
  double r = x->data[0];
  double s = x->data[1];
  double t = x->data[2];
  if((r<0) || (s<0) || (t<0) )
    return log(0.0);
  return loglik(x)
    + log(gsl_ran_gamma_pdf(r, 1.0, 1.0))
    + log(gsl_ran_gamma_pdf(s, 1.0, 1.0))
    + log(gsl_ran_gamma_pdf(t, 1.0, 1.0));
}

int main(int argc, char** argv) {
  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  /*set observed data*/
  y_obs = gsl_vector_alloc(S);
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
      gsl_vector_set(y_obs, i*3 + j, gsl_ran_gaussian(r, 1.0));
  }

  /*read spatial info*/
  gsl_matrix* XY = gsl_matrix_alloc(S, 2);
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
      gsl_matrix_set(XY, i*3 + j, 0, (double) i);
      gsl_matrix_set(XY, i*3 + j, 1, (double) j);
    }
  gsl_matrix* D = gsl_matrix_alloc(S, S);
  mcmclib_spatial_distances(D, XY);
  gsl_matrix_free(XY);
  gsl_vector* mu = gsl_vector_alloc(S);
  gsl_vector_set_all(mu, 0.0);
  theta = gsl_vector_alloc(3);
  rho_v = gsl_vector_subvector(theta, 0, 1);
  sigma_v = gsl_vector_subvector(theta, 1, 1);
  tausq_v = gsl_vector_subvector(theta, 2, 1);
  lik = mcmclib_spatial_lpdf_alloc(mu, &(rho_v.vector), &(sigma_v.vector), &(tausq_v.vector), D);

  /*set starting chain values*/
  int d = theta->size;
  gsl_vector* x = gsl_vector_alloc(d);
  gsl_vector_set_all(x, 0.5);

  /*alloc a new sampler*/
  mcmclib_mh* sampler = mcmclib_gauss_rw_alloc(r, target_logdensity, NULL, x, V0);

  /*open output file*/
  FILE* out = fopen(OUTPUT_FILE, "w");
  /*print out csv header*/
  fprintf(out, "rho, sigma, tausq\n");

  /*main MCMC loop*/
  for(int i=0; i<N; i++) {
    mcmclib_mh_update(sampler);
    for(int j=0; j<(d-1); j++)
      fprintf(out, "%f, ", gsl_vector_get(x, j));
    fprintf(out, "%f\n", gsl_vector_get(x, d-1));
  }

  fclose(out);
  gsl_rng_free(r);
  gsl_vector_free(x);
  mcmclib_gauss_rw_free(sampler);
  mcmclib_spatial_lpdf_free(lik);
  gsl_matrix_free(D);
  gsl_vector_free(mu);
  gsl_vector_free(theta);
  gsl_vector_free(y_obs);
}
