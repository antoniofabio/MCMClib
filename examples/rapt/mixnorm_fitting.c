#include <gsl/gsl_math.h>
#include "mixnorm_fitting.h"

mcmclib_mixolem_suff* gamma_hat;
FILE *out_beta, *out_mu, *out_Sigma;

/*init fitted mixture data*/
void mixnorm_fitting_init(gsl_rng* r) {
  gamma_hat = mcmclib_mixolem_suff_alloc(K, DIM);
  gsl_vector_set_all(gamma_hat->delta, 1.0 / (double) K);
  for(int k=0; k<K; k++) {
    gsl_vector_set_all(gamma_hat->delta_x[k], MU0 * pow(-1.0, k+1) * (0.75 + 0.5 * gsl_rng_uniform(r)));
    gsl_matrix_set_identity(gamma_hat->delta_xx[k]);
    gsl_matrix_scale(gamma_hat->delta_xx[k], V0[k] * (0.75 + 0.5 * gsl_rng_uniform(r)));
  }

  out_beta = fopen("beta_hat.dat", "w");
  out_mu = fopen("mu_hat.dat", "w");
  out_Sigma = fopen("Sigma_hat.dat", "w");
}

void mixnorm_fitting_free() {
  mcmclib_mixolem_suff_free(gamma_hat);
  fclose(out_beta);
  fclose(out_mu);
  fclose(out_Sigma);
}

/*save currently estimated parameters*/
void mixnorm_fitting_store() {
  gsl_vector_fprintf(out_beta, gamma_hat->delta, "%f");
  for(int k=0; k<K; k++){
    gsl_vector_fprintf(out_mu, gamma_hat->delta_x[k], "%f");
    gsl_matrix_fprintf(out_Sigma, gamma_hat->delta_xx[k], "%f");
  }
}
