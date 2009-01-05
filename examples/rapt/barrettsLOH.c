/**
OLEM-RAPT, Barretts LOH data
*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
/*recursive mixture fitting*/
#include <mixem_online.h>
/*OLEM-RAPT algorithm*/
#include <olemrapt_inca.h>

/** input parameters*/
/*total number of iterations*/
const int N = 50000;
/*number of parallel chains to run*/
const int K = 20;
/*HCL burn in*/
const int T0 = 1000;
/*starting means guesses*/
double M0[] = {-5.0, -5.0, 5.0, 5.0,
		     5.0, 5.0, -5.0, -5.0};
/*starting variance guesses*/
double V0[] = {4.0, 4.0};
/*starting overall variance*/
const double S0 = 6.0;
/*target space dimension*/
const int DIM = 4;
/*scaling factor*/
#define SCALING_FACTOR ((2.38 * 2.38) / (double) DIM)
/*data filenamee*/
#define DATA_FNAME "BarrettsLOH.dat"
/*no. of rows*/
const int NR = 40;

gsl_matrix* LOH;

static double ldbb2(int x, int size, double prob, double omega) {
  if(omega < 1e-5)
    return(log(gsl_ran_binomial_pdf(x, prob, size)));

  double theta = 1.0 / omega;
  double ans = gsl_sf_lnchoose(size, x) - gsl_sf_lnbeta(theta * (1-prob), theta * prob) +
    gsl_sf_lnbeta(size - x + theta * (1-prob), x + theta * prob);
  return(ans);
}

static double my_lf(int x, int n, double eta, double pi1, double pi2, double gamma) {
  double omega2 = exp(gamma) / (2.0 * (1.0 + exp(gamma)));
  double ans = eta * gsl_ran_binomial_pdf(x, pi1, n) +
    (1.0 - eta) * exp(ldbb2(x, n, pi2, omega2));
  return log(ans);
}

static void log_pi(double* theta, double* ans) {
  double eta = theta[0];
  double pi1 = theta[1];
  double pi2 = theta[2];
  double gamma = theta[3];
  if((eta < 0.01)  || (pi1 < 0.01) || (pi2 < 0.01) ||
     (eta > 0.99) || (pi1 > 0.99) || (pi2 > 0.99) ||
     (abs(gamma) >= 30.0)) {
    ans[0] = log(0.0);
    return;
  }
  ans[0] = 0.0;
  for(int i=0; i<NR; i++)
    ans[0] += my_lf(gsl_matrix_get(LOH, i, 0), gsl_matrix_get(LOH, i, 1), theta[0], theta[1], theta[2], theta[3]);
}

/*target distribution*/
double target_logdensity(void* ignore, gsl_vector* x) {
  double ans;
  log_pi(x->data, &ans);
  return ans;
}

#define STORE_VALUES \
  for(int i=0; i<K; i++) { \
    for(int j=0; j<(d-1); j++) \
      fprintf(out_X, "%f, ", gsl_vector_get(xx[i], j)); \
    fprintf(out_X, "%f\n", gsl_vector_get(xx[i], d-1)); \
  }

int main(int argc, char** argv) {
  /*load and setup input data*/
  LOH = gsl_matrix_alloc(NR, 2);
  FILE* fdata = fopen(DATA_FNAME, "r");
  if(!fdata)
    return(1);
  gsl_matrix_fscanf(fdata, LOH);
  fclose(fdata);

  int d = DIM;
  /*set starting mixture parameters values*/
  gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
  gsl_matrix_set_identity(Sigma_zero);
  gsl_matrix_scale(Sigma_zero, S0);
  /*means*/
  gsl_vector_view m0v[2];
  m0v[0] = gsl_vector_view_array(M0, DIM);
  m0v[1] = gsl_vector_view_array(M0 + DIM, DIM);
  gsl_vector* mu0[2];
  mu0[0] = &(m0v[0].vector);
  mu0[1] = &(m0v[1].vector);
  /*variances*/
  gsl_matrix* Sigma[2];
  for(int k=0; k<2; k++) {
    Sigma[k] = gsl_matrix_alloc(d, d);
    gsl_matrix_set_identity(Sigma[k]);
    gsl_matrix_scale(Sigma[k], V0[k]);
  }
  /*weights*/
  gsl_vector* beta = gsl_vector_alloc(2);
  gsl_vector_set_all(beta, 0.5);

  /*alloc a new RNG*/
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* xx[K];
  /*set starting values at random*/
  for(int k=0; k<K; k++) {
    xx[k] = gsl_vector_alloc(d);
    do {
      for(int j=0; j<(d-1); j++)
	gsl_vector_set(xx[k], j, gsl_rng_uniform(r));
      gsl_vector_set(xx[k], d-1, gsl_rng_uniform(r) * 60.0 - 30.0);
    } while(!isfinite(target_logdensity(NULL, xx[k])));
  }
  /*alloc a new OLEM-RAPT INCA sampler*/
  mcmclib_olemrapt_inca* sampler =
    mcmclib_olemrapt_inca_alloc(r, target_logdensity, NULL,
				xx, T0, Sigma_zero,
				beta, mu0, Sigma, K);

  /*open output files*/
  FILE* out_X = fopen("barrettsLOH_X.csv", "w");
  /*print out csv header*/
  fprintf(out_X, "eta, pi1, pi2, gamma\n");
  STORE_VALUES;
  for(int i=0; i<K; i++)
    printf("lpi(%d) = %f\n", i, target_logdensity(NULL, xx[i]));

  /*main MCMC loop*/
  for(int n=0; n<N; n++) {
    /*update chain value*/
    mcmclib_olemrapt_inca_update(sampler);
    mcmclib_olemrapt_inca_update_proposals(sampler);
    /*store sampled values*/
    STORE_VALUES;
  }

  /*free resources*/
  fclose(out_X);
  gsl_matrix_free(Sigma_zero);
  gsl_rng_free(r);
  mcmclib_olemrapt_inca_free(sampler);
  for(int m=0; m<K; m++)
    gsl_vector_free(xx[m]);

  return 0;
}
