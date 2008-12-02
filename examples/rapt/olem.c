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
/*target distribution*/
#include "mixnorm_target.h"
#include "mixnorm_fitting.h"

/** input parameters*/
/**number of parallel chains*/
const int M = 4;
/*chain length*/
const int N = 25000;
/*burn in length*/
const int T0 = 200;
/*scaling factor*/
#define SCALING_FACTOR ((2.38 * 2.38) / (double) DIM)
/*starting global variance guess*/
#define SIGMA0_GLOBAL (pow(MU0 * 2.0, 2.0) + V0[0] + V0[1])
/***/

int main(int argc, char** argv) {
  /*init target distribution data*/
  mcmclib_mixnorm_lpdf* pi_target = mixnorm_target_alloc();

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  /*init fitted mixture data*/
  mixnorm_fitting_init(r);

  /*set starting guess metropolis covariance matrices*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, SIGMA0_GLOBAL * SCALING_FACTOR);

  /*open output file*/
  FILE* out_X = fopen("X.dat", "w"); /*sampled data*/

  /*init current chains values*/
  gsl_vector* x[M];
  for(int m=0; m<M; m++) {
    x[m] = gsl_vector_alloc(DIM);
    gsl_vector_set_all(x[m], MU0 * pow(-1.0, m));
  }

  /*alloc a new set of OLEM-RAPT INCA samplers*/
  mcmclib_olemrapt_inca* sampler =
    mcmclib_olemrapt_inca_alloc(r,
				mcmclib_mixnorm_lpdf_compute, pi_target,
				x,
				T0, Sigma_zero,
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
    /*store estimated mixture parameters values*/
    mixnorm_fitting_store();
  }

  /*free resources*/
  fclose(out_X);
  mixnorm_target_free(pi_target);
  gsl_rng_free(r);
  mcmclib_olemrapt_inca_free(sampler);
  gsl_matrix_free(Sigma_zero);
  for(int m=0; m<M; m++)
    gsl_vector_free(x[m]);
  mixnorm_fitting_free();
}
