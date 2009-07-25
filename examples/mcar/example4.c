/**Adaptive Gaussian Random Walk and RAPTOR on a Poisson model
   with MCAR effects example. Real data. */
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <mvnorm.h>
#include <mcar_tilde.h>
#include <mcar_model.h>
#include <pois_model.h>
#include <gauss_am.h>
#include <raptor.h>

/* P=3, DIM=95: 0.06590 secs per iteration */
#define N 50000
#define THIN 10
#define T0 15000
#define V0 0.4
#define SF 0.05 /*scaling factor*/

#define P 3
#define DIM 95

gsl_rng* rng;
gsl_vector *alpha12sigma, *alphasigmag;
const int a12s_size = P * P;
const int as_size = P * (P-1) / 2 + P;
mcmclib_mcar_tilde_lpdf* mcar_lpdf;
mcmclib_mcar_model* mcar_model;
gsl_vector* mcar_phi;
gsl_vector* denom;
gsl_vector* offset;
gsl_vector* y;
gsl_matrix* X;
mcmclib_pmodel_sampler* model;
mcmclib_amh* sampler[3 + DIM];
int util_j[DIM];
gsl_vector_view phij[DIM];

/*mixture parameters*/
gsl_vector *a12s_beta, *as_beta;
gsl_vector *a12s_mu[2], *as_mu[2];
gsl_matrix *a12s_Sigma[2], *as_Sigma[2];

gsl_vector_view asmu_view[2];
double asmu_data[2][6] = {{0.6402105,  1.6130955, -0.1238163, -0.1596428, -0.9227976, -2.6514444},
			  {-0.4479953, 0.7868348, -0.3515181, -0.3648961, -1.0954467, -1.8379377}};

gsl_vector_view a12smu_view[2];
double a12smu_data[2][9] = {{-0.2546310, -0.3525798, 0.1891869, -0.2442633, 0.3454967, -0.0238164,
			     -0.1346993, -0.3288677, -1.5137814},
			    {-0.28717678, -0.18210887,  0.07761313, -0.60640886, -0.07559915, -0.11495301,
			     -0.42148466, -1.24859913, -2.29223423}};

static void init_mixture_params() {
  a12s_beta = gsl_vector_alloc(2);
  gsl_vector_set_all(a12s_beta, 0.5);
  as_beta = gsl_vector_alloc(2);
  gsl_vector_set_all(as_beta, 0.5);
  for(int k=0; k<2; k++) {
    a12smu_view[k] = gsl_vector_view_array(a12smu_data[k], a12s_size);
    a12s_mu[k] = &(a12smu_view[k].vector);
    a12s_Sigma[k] = gsl_matrix_alloc(a12s_size, a12s_size);
    gsl_matrix_set_identity(a12s_Sigma[k]);
    asmu_view[k] = gsl_vector_view_array(asmu_data[k], as_size);
    as_mu[k] = &(asmu_view[k].vector);
    gsl_vector_set_all(as_mu[k], 2.0);
    as_Sigma[k] = gsl_matrix_alloc(as_size, as_size);
    gsl_matrix_set_identity(as_Sigma[k]);
  }
}

static void free_mixture_params() {
  gsl_vector_free(a12s_beta);
  gsl_vector_free(as_beta);
  for(int k=0; k<2; k++) {
    gsl_matrix_free(a12s_Sigma[k]);
    gsl_matrix_free(as_Sigma[k]);
  }
}

static void update_offset() {
  gsl_vector_memcpy(offset, denom);
  gsl_vector_add(offset, mcar_phi);
}

static double phij_lpdf(void* pj, gsl_vector* phij) {
  int j = ((int*) pj)[0];
  double prior = mcmclib_mcar_model_phi_fcond(mcar_model, j, phij);
  gsl_vector* phi_old = gsl_vector_alloc(DIM * P);
  gsl_vector_memcpy(phi_old, mcar_phi);
  gsl_vector_view phiv = gsl_vector_subvector(mcar_phi, j*P, P);
  gsl_vector_memcpy(&phiv.vector, phij);
  update_offset();
  double lik = mcmclib_pois_model_llik(model->model, model->sampler->mh->x);
  gsl_vector_memcpy(mcar_phi, phi_old);
  update_offset();
  gsl_vector_free(phi_old);
  return prior + lik;
}

void init_chains() {
  alpha12sigma = mcar_lpdf->alpha12sigma;
  gsl_vector_set_all(alpha12sigma, -1.0);
  gsl_matrix* Sigma0 = gsl_matrix_alloc(a12s_size, a12s_size);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / ((double)a12s_size));
  sampler[0] = mcmclib_raptor_alloc(rng, mcmclib_mcar_model_alpha12sigma_lpdf,
				    mcar_model, alpha12sigma, T0,
				    Sigma0, a12s_beta, a12s_mu, a12s_Sigma);
  mcmclib_raptor_set_sf_local(sampler[0], SF * 1.5);
  mcmclib_raptor_set_sf_global(sampler[0], SF);
  mcmclib_raptor_set_alpha(sampler[0], 0.8);
  gsl_matrix_free(Sigma0);

  alphasigmag = mcar_lpdf->alphasigmag;
  gsl_vector_set_all(alphasigmag, -1.0);
  Sigma0 = gsl_matrix_alloc(as_size, as_size);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / (double) as_size);
  sampler[1] = mcmclib_raptor_alloc(rng, mcmclib_mcar_model_alphasigma_lpdf,
				    mcar_model, alphasigmag, T0,
				    Sigma0, as_beta, as_mu, as_Sigma);
  mcmclib_raptor_set_sf_local(sampler[1], SF);
  mcmclib_raptor_set_sf_global(sampler[1], SF);
  mcmclib_raptor_set_alpha(sampler[1], 0.8);
  gsl_matrix_free(Sigma0);

  model = mcmclib_pmodel_sampler_alloc(X, y, offset, rng, 1e-3, T0);
  gsl_vector_set_all(model->model->beta, 0.0);
  sampler[2] = model->sampler;

  Sigma0 = gsl_matrix_alloc(P, P);
  gsl_matrix_set_identity(Sigma0);
  gsl_matrix_scale(Sigma0, V0 / (double) (P));
  for(int j=0; j<DIM; j++) {
    util_j[j] = j;
    phij[j] = gsl_vector_subvector(mcar_phi, j*P, P);
    sampler[j+3] = mcmclib_gauss_am_alloc(rng, phij_lpdf, util_j + j, &phij[j].vector,
					  Sigma0, T0);
  }
  gsl_matrix_free(Sigma0);
}

void free_chains() {
  for(int i=0; i<2; i++)
    mcmclib_raptor_free(sampler[i]);
  mcmclib_pmodel_sampler_free(model);
  gsl_rng_free(rng);
}

int main(int argc, char** argv) {
  rng = gsl_rng_alloc(gsl_rng_default);

  /* model setup */
  gsl_matrix* W = gsl_matrix_alloc(DIM, DIM);
  FILE* in_W = fopen("W_2.dat", "r");
  gsl_matrix_fscanf(in_W, W);
  fclose(in_W);
  mcar_lpdf = mcmclib_mcar_tilde_lpdf_alloc(P, W);
  gsl_matrix_free(W);
  mcar_phi = gsl_vector_alloc(P * DIM);
  gsl_vector_set_zero(mcar_phi);
  denom = gsl_vector_alloc(P * DIM);
  offset = gsl_vector_alloc(P * DIM);
  FILE* in_offset = fopen("offset_2.dat", "r");
  gsl_vector_fscanf(in_offset, denom);
  fclose(in_offset);
  for(int i=0; i < P*DIM; i++)
    gsl_vector_set(denom, i, log(gsl_vector_get(denom, i)));
  update_offset();
  mcar_model = mcmclib_mcar_model_alloc(mcar_lpdf, mcar_phi);

  y = gsl_vector_alloc(P * DIM);
  FILE* in_y = fopen("y_2.dat", "r");
  gsl_vector_fscanf(in_y, y);
  fclose(in_y);
  X = gsl_matrix_alloc(P * DIM, P);
  gsl_matrix_set_zero(X);
  for(int d=0; d<DIM; d++)
    for(int i=0; i<P; i++)
      gsl_matrix_set(X, d*P + i, i, 1.0);

  init_mixture_params();
  init_chains();
  FILE* out_a12s = fopen("chain_alpha12sigma.dat", "w");
  FILE* out_as = fopen("chain_alphasigma.dat", "w");
  FILE* out_beta = fopen("chain_beta.dat", "w");
  FILE* out_lpdf = fopen("chain_lpdf.dat", "w");
  FILE* out_gii = fopen("chain_gammaii.dat", "w");
  FILE* out_bii = fopen("chain_bii.dat", "w");
  gsl_vector_view gii_v = gsl_matrix_diagonal(mcar_lpdf->Gamma);
  gsl_vector_view bii_v = gsl_matrix_diagonal(mcar_lpdf->B_tilde);
  for(int i=0; i<N; i++) {
    if (( (i+1) % THIN ) == 0) {
      gsl_vector_fprintf(out_a12s, alpha12sigma, "%f");
      gsl_vector_fprintf(out_as, alphasigmag, "%f");
      gsl_vector_fprintf(out_beta, sampler[2]->mh->x, "%f");
      fprintf(out_lpdf, "%f\n",
	      mcmclib_mcar_model_alphasigma_lpdf(mcar_model, alphasigmag));
      gsl_vector_fprintf(out_gii, &gii_v.vector, "%f");
      gsl_vector_fprintf(out_bii, &bii_v.vector, "%f");
      fflush(out_beta);
      fflush(out_lpdf);
      fflush(out_a12s);
      fflush(out_as);
      fflush(out_gii);
      fflush(out_bii);
    }
    for(int j=0; j<2; j++) {
      mcmclib_amh_update(sampler[j]);
      assert(gsl_finite(sampler[j]->mh->logdistr_old));
      mcmclib_mcar_tilde_lpdf_update_vcov(mcar_lpdf);
    }
    for(int j=2; j<DIM+3; j++) {
      mcmclib_amh_update(sampler[j]);
      assert(gsl_finite(sampler[j]->mh->logdistr_old));
      update_offset();
    }
  }
  fclose(out_a12s);
  fclose(out_as);
  fclose(out_beta);
  fclose(out_lpdf);
  fclose(out_gii);
  fclose(out_bii);

  //  free_chains();
  free_mixture_params();
  mcmclib_mcar_model_free(mcar_model);
  gsl_vector_free(mcar_phi);
  mcmclib_mcar_tilde_lpdf_free(mcar_lpdf);
}