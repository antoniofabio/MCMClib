/***********/
/*OLEM-RAPT*/
/***********/
#include <stdio.h>
#include <memory.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/*multivariate normal and mixture of normals distributions*/
#include <mvnorm.h>
#include <mixnorm.h>
/*OLEM-RAPT algorithm*/
#include <olemrapt_inca.h>

#include "mixnorm_target.h"

/** input parameters*/
/**number of parallel chains*/
#define M 4
/*chain length*/
static int N = 25000;
/*burn in length*/
#define T0 ((int) pow(DIM, 2.0) * 25)
/*scaling factor*/
#define SCALING_FACTOR (2.38 * 2.38 / (double) DIM)
/*starting global variance guess*/
#define SIGMA0_GLOBAL (pow(MU0 * 2.0, 2.0) + V0[0] + V0[1])
/***/

/********************************/
/*boundary function support data*/
/********************************/
/*estimated mixture par. values*/
mcmclib_mixolem_suff* gamma_hat;
/**/

/*output files handlers*/
FILE *out_X, *out_beta, *out_mu, *out_Sigma;

void save_gamma_hat();
int main(int argc, char** argv) {
  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  /*current chains values*/
  gsl_vector* x[M];
  for(int m=0; m<M; m++) {
    x[m] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(x[m], MU0 * pow(-1.0, m));
  }
  /*init target distribution data*/
  mcmclib_mixnorm_lpdf* pi_target = mixnorm_target_alloc();

  /*init fitted mixture data*/
  gamma_hat = mcmclib_mixolem_suff_alloc(K, DIM);
  gsl_vector_set_all(gamma_hat->delta, 1.0 / (double) K);
  for(int k=0; k<K; k++) {
    gsl_vector_set_all(gamma_hat->delta_x[k], MU0 * pow(-1.0, k+1) * (0.75 + 0.5 * gsl_rng_uniform(r)));
    gsl_matrix_set_identity(gamma_hat->delta_xx[k]);
    gsl_matrix_scale(gamma_hat->delta_xx[k], V0[k] * (0.75 + 0.5 * gsl_rng_uniform(r)));
  }

  /*set starting guess metropolis covariance matrices*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, SIGMA0_GLOBAL * SCALING_FACTOR);

  /*open output files*/
  out_X = fopen("X.dat", "w"); /*sampled data*/
  out_beta = fopen("beta_hat.dat", "w"); /*beta_hat*/
  out_mu = fopen("mu_hat.dat", "w"); /*mu_hat*/
  out_Sigma = fopen("Sigma_hat.dat", "w"); /*Sigma_hat*/

  /*alloc a new OLEM-RAPT INCA samplers*/
  mcmclib_olemrapt_inca* sampler =
    mcmclib_olemrapt_inca_alloc(r,
				mcmclib_mixnorm_lpdf_compute, pi_target,
				x, T0, Sigma_zero,
				gamma_hat->delta,
				gamma_hat->delta_x,
				gamma_hat->delta_xx,
				M);

  /*main MCMC loop*/
  for(int n=0; n<N; n++) {
    /*update chain value*/
    mcmclib_olemrapt_inca_update(sampler);
    mcmclib_olemrapt_inca_update_proposals(sampler);
    /*store sampled values*/
    for(int m=0; m<M; m++)
      gsl_vector_fprintf(out_X, x[m], "%f");
    save_gamma_hat();
  }

  fclose(out_X);
  fclose(out_mu);
  fclose(out_Sigma);
  fclose(out_beta);
  mixnorm_target_free(pi_target);
  gsl_rng_free(r);
  mcmclib_olemrapt_inca_free(sampler);
  gsl_matrix_free(Sigma_zero);
}

/*save currently estimated parameters*/
void save_gamma_hat() {
  gsl_vector_fprintf(out_beta, gamma_hat->delta, "%f");
  for(int k=0; k<K; k++){
    gsl_vector_fprintf(out_mu, gamma_hat->delta_x[k], "%f");
    gsl_matrix_fprintf(out_Sigma, gamma_hat->delta_xx[k], "%f");
  }
}
