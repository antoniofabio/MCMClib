#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gauss_am.h>

#define OUTPUT_FILE "data.csv"
#define N 100000
#define T0 20
#define V0 0.1

/*target distribution: uniform in the unit cube*/
double target_logdensity(void* ignore, gsl_vector* x) {
	int d = x->size;
	for(int i=0; i<d; i++) {
		double xi = gsl_vector_get(x, i);
		if((xi < 0.0) || (xi > 1.0))
			return log(0.0);
	}
	return 0;
}

int main(int argc, char** argv) {
	int d = 3;
	/*set starting guess covariance matrix*/
	gsl_matrix* Sigma_zero = gsl_matrix_alloc(d, d);
	gsl_matrix_set_identity(Sigma_zero);
	gsl_matrix_scale(Sigma_zero, V0);
	
	/*alloc a new adaptive sampler*/
	mcmclib_gauss_am* sampler = mcmclib_gauss_am_alloc(Sigma_zero, T0);
	/*alloc a new RNG*/
	gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
	/*alloc current chain value space*/
	gsl_vector* x = gsl_vector_alloc(d);
	/*set starting value*/
	gsl_vector_set_all(x, 0.5);

	/*open output file*/
	FILE* out = fopen(OUTPUT_FILE, "w");
	/*print out csv header*/
	fprintf(out, "x1, x2, x3\n");

	/*main MCMC loop*/
	for(int i=0; i<N; i++) {
		mcmclib_gauss_am_update(sampler, r, target_logdensity, x, NULL);
		for(int j=0; j<(d-1); j++)
			fprintf(out, "%f, ", gsl_vector_get(x, j));
		fprintf(out, "%f\n", gsl_vector_get(x, d-1));
	}

	fclose(out);
	gsl_rng_free(r);
	mcmclib_gauss_am_free(sampler);
	gsl_matrix_free(Sigma_zero);
}
