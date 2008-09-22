/**
INCA example, Barretts LOH data
*/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gauss_inca.h>

/*total number of iterations*/
#define N 50000
/*number of parallel chains to run*/
#define K 5
/*HCL burn in*/
#define T0 100
/*starting variance guess*/
#define V0 1.0
/*target space dimension*/
#define DIM 4
/*data filenamee*/
#define DATA_FNAME "BarrettsLOH.dat"
/*no. of rows*/
#define NR 40

gsl_vector* xx;
gsl_vector* nn;

/*no. of combinations of 'n' over 'r' (log)*/
double lnchoose(int n, int r){
	double num = 0.0;
	double den = 0.0;
	while(r) {
		num += log(n--);
		den += log(r--);
	}
	return num-den;
}

double my_lngamma(double x) {
	gsl_sf_result result;
	int status = gsl_sf_lngamma_e (x, &result);
	if(status == GSL_SUCCESS)
		return(result.val);
	else
		return(1.0 / 0.0);
}

double f(int in_x, int in_n, double eta, double pi1, double pi2, double gamma) {
	double omega2 = exp(gamma) / (2.0 * (1.0 + exp(gamma)));
	double a = 0.0;
	double b = 0.0;
	double x = in_x;
	double n = in_n;
	a = lnchoose(n, x) + x * log(pi1) + (n-x) * log(1.0 - pi1);
	b = lnchoose(n, x) + my_lngamma(1.0 / omega2) + my_lngamma(x + pi2 / omega2);
	b -= my_lngamma(pi2 / omega2) + my_lngamma((1.0 - pi2) / omega2) + my_lngamma(n - x + (1.0 - pi2) / omega2) + my_lngamma(n + 1.0 / omega2);
	double ans = eta * exp(a) + (1.0 - eta) * exp(b);
	return log(ans);
}

/*target distribution*/
double target_logdensity(void* ignore, gsl_vector* x) {
	double eta = gsl_vector_get(x, 0);
	double pi1 = gsl_vector_get(x, 1);
	double pi2 = gsl_vector_get(x, 2);
	double gamma = gsl_vector_get(x, 3);
	if((eta < 0.0)  || (pi1 < 0.0) || (pi2 < 0.0) ||
			(eta > 1.0) || (pi1 > 1.0) || (pi2 > 1.0) ||
			(abs(gamma) >= 30.0)) {
		return log(0.0);
	}
	double ans = 0.0;
	for(int i=0; i<NR; i++)
		ans += f(gsl_vector_get(xx,i), gsl_vector_get(nn,i), eta, pi1, pi2, gamma);
	return ans;
}

int main(int argc, char** argv) {
	gsl_set_error_handler_off ();

	/*load and setup input data*/
	gsl_matrix* data = gsl_matrix_alloc(NR, 2);
	FILE* fdata = fopen(DATA_FNAME, "r");
	if(!fdata)
		return(1);
	gsl_matrix_fscanf(fdata, data);
	fclose(fdata);
	gsl_vector_view xv = gsl_matrix_column(data, 0);
	xx = &(xv.vector);
	gsl_vector_view nv = gsl_matrix_column(data, 1);
	nn = &(nv.vector);

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
