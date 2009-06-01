/**Adaptive Gaussian Random Walk on an MCAR model example*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>
#include <mcar_tilde.h>
#include <mcar_model.h>
#include <matrix.h>

/*P=3, DIM=95: 0.07379 secs per iteration
  P=6, DIM=10: 0.00232 secs per iteration
  P=3, DIM=10: 0.00039 secs per iteration
  P=3, DIM=5:  0.00016 secs per iteration
*/
#define N 10000
#define THIN 10
#define T0 5000
#define V0 0.4
#define SF 0.2

#define P 3
#define DIM 95

gsl_vector *alpha12sigma, *alphasigmag;
gsl_rng* rng;
mcmclib_mcar_tilde_lpdf* lpdf;
mcmclib_mcar_model* model;
mcmclib_amh* sampler[2];

void init_chains() {
  alpha12sigma = lpdf->alpha12sigma;
  gsl_vector_set_all(alpha12sigma, -1.0);
  gsl_matrix* Sigma0 = gsl_matrix_alloc(P*P, P*P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / ((double)(P * P)));
  sampler[0] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alpha12sigma_lpdf,
				      model, alpha12sigma, Sigma0, T0);
  ((mcmclib_gauss_am_suff*) sampler[0]->suff)->sf = SF;
  gsl_matrix_free(Sigma0);

  alphasigmag = lpdf->alphasigmag;
  gsl_vector_set_all(alphasigmag, -1.0);
  Sigma0 = gsl_matrix_alloc(P*(P-1)/2 + P, P*(P-1)/2 + P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / (double) (P*(P-1)/2 + P));
  sampler[1] = mcmclib_gauss_am_alloc(rng, mcmclib_mcar_model_alphasigma_lpdf,
				      model, alphasigmag, Sigma0, T0);
  ((mcmclib_gauss_am_suff*) sampler[1]->suff)->sf = SF;
  gsl_matrix_free(Sigma0);
}

void free_chains() {
  for(int i=0; i<2; i++)
    mcmclib_gauss_am_free(sampler[i]);
  gsl_rng_free(rng);
}

int main(int argc, char** argv) {
  rng = gsl_rng_alloc(gsl_rng_default);

  /* model setup */
  gsl_matrix* W = gsl_matrix_alloc(DIM, DIM);
  gsl_matrix_set_zero(W);
  for(int i=0; i<(DIM-1); i++) {
    gsl_matrix_set(W, i, i+1, 1.0);
    gsl_matrix_set(W, i+1, i, 1.0);
  }
  lpdf = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  gsl_matrix_free(W);
  gsl_vector* e = gsl_vector_alloc(P * DIM);
  FILE* phi_in = fopen("phi.dat", "r");
  gsl_vector_fscanf(phi_in, e);
  fclose(phi_in);
  model = mcmclib_mcar_model_alloc(lpdf, e);

  init_chains();
  FILE* out_a12s = fopen("chain_alpha12sigma.dat", "w");
  FILE* out_as = fopen("chain_alphasigma.dat", "w");
  FILE* out_lpdf = fopen("chain_lpdf.dat", "w");
  for(int i=0; i<N; i++) {
    if (( (i+1) % THIN ) == 0) {
      gsl_vector_fprintf(out_a12s, alpha12sigma, "%f");
      gsl_vector_fprintf(out_as, alphasigmag, "%f");
      fprintf(out_lpdf, "%f\n", mcmclib_mcar_model_alpha12sigma_lpdf(model, alpha12sigma));
    }
    for(int j=0; j<2; j++) {
      mcmclib_amh_update(sampler[j]);
      assert(gsl_finite(sampler[j]->mh->logdistr_old));
    }
  }
  fclose(out_a12s);
  fclose(out_as);
  fclose(out_lpdf);

  free_chains();
  mcmclib_mcar_model_free(model);
  gsl_vector_free(e);
  mcmclib_mcar_tilde_lpdf_free(lpdf);
}
