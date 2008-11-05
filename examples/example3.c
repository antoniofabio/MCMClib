/**
INCA example, Barretts LOH data
*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gauss_inca.h>

/*total number of iterations*/
#define N 50000
/*number of parallel chains to run*/
#define K 20
/*HCL burn in*/
#define T0 1000
/*starting variance guess*/
#define V0 0.5
/*target space dimension*/
#define DIM 4
/*data filenamee*/
#define DATA_FNAME "BarrettsLOH.dat"
/*no. of rows*/
#define NR 40

gsl_matrix* LOH;
gsl_vector* LOH_x;
gsl_vector* LOH_n;

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

int main(int argc, char** argv) {
	/*load and setup input data*/
	LOH = gsl_matrix_alloc(NR, 2);
	FILE* fdata = fopen(DATA_FNAME, "r");
	if(!fdata)
		return(1);
	gsl_matrix_fscanf(fdata, LOH);
	fclose(fdata);
/*	gsl_vector_view xv = gsl_matrix_column(data, 0);
	LOH_x = &(xv.vector);
	gsl_vector_view nv = gsl_matrix_column(data, 1);
	LOH_n = &(nv.vector);*/

	int d = DIM;
	/*set starting guess covariance matrix*/
	gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
	gsl_matrix_set_identity(Sigma_zero);
	gsl_matrix_scale(Sigma_zero, V0);

	/*alloc a new RNG*/
	gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
	/*lets go sample from K different chains in parallel*/
	mcmclib_gauss_inca_pool* pool = mcmclib_gauss_inca_pool_alloc(Sigma_zero, T0, K);
	gsl_vector* xx[K];
	mcmclib_gauss_inca* sampler[K];
	for(int k=0; k<K; k++) {
		sampler[k] = mcmclib_gauss_inca_alloc(pool);
		xx[k] = gsl_vector_alloc(d);
		/*set starting values at random*/
		do {
			for(int j=0; j<(d-1); j++)
				gsl_vector_set(xx[k], j, gsl_rng_uniform(r));
			gsl_vector_set(xx[k], d-1, gsl_rng_uniform(r)*60.0 - 30.0);
		} while(!isfinite(target_logdensity(NULL, xx[k])));
	}

	/*open output files*/
	FILE* out[K];
	char filename[512];
	for(int i=0; i<K; i++) {
		sprintf(filename, "chain_%d.csv", i);
		out[i] = fopen(filename, "w");
		/*print out csv header*/
		fprintf(out[i], "eta, pi1, pi2, gamma\n");
		for(int j=0; j<(d-1); j++)
			fprintf(out[i], "%f, ", gsl_vector_get(xx[i], j));
		fprintf(out[i], "%f\n", gsl_vector_get(xx[i], d-1));

		printf("lpi(%d) = %f\n", i, target_logdensity(NULL, xx[i]));
	}

	/*main MCMC loop: for each time step iterate troughout the K chains*/
	for(int i=0; i<N; i++) for(int k=0; k<K; k++) {
		mcmclib_gauss_inca_update(sampler[k], r, target_logdensity, xx[k], NULL);
		for(int j=0; j<(d-1); j++)
			fprintf(out[k], "%f, ", gsl_vector_get(xx[k], j));
		fprintf(out[k], "%f\n", gsl_vector_get(xx[k], d-1));
	}

	/*release system resources*/
	for(int k=0; k<K; k++) {
		fclose(out[k]);
		mcmclib_gauss_inca_free(sampler[k]);
		gsl_vector_free(xx[k]);
	}
	mcmclib_gauss_inca_pool_free(pool);
	gsl_rng_free(r);
	gsl_matrix_free(Sigma_zero);

	return 0;
}
