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
#define N 100000
/*number of parallel chains to run*/
#define K 3
/*HCL burn in*/
#define T0 100
/*starting variance guess*/
#define V0 0.1
/*target space dimension*/
#define DIM 4
/*data filenamee*/
#define DATA_FNAME "BarrettsLOH.dat"
/*no. of rows*/
#define NR 40

gsl_vector* xx;
gsl_vector* nn;

double choose(int n, int r){
	double num = 1.0;
	double den = 1.0;
	while(r-1) {
		num *= n--;
		den *= r--;
	}
	return num/den;
}

double f(int x, int n, double eta, double pi1, double pi2, double gamma) {
	eta = gsl_cdf_logistic_P(eta, 1.0);
	pi1 = gsl_cdf_logistic_P(pi1, 1.0);
	pi2 = gsl_cdf_logistic_P(pi2, 1.0);
	double omega2 = exp(gamma) / (2 * (1 + exp(gamma)));
	double a = 0.0;
	double b = 0.0;
	a = choose(n, x) * pow(pi1, x) * pow(1 - pi1, n - x);
	b = choose(n, x) * gsl_sf_gamma(1/omega2) * gsl_sf_gamma(x + pi2 / omega2);
	b /= gsl_sf_gamma(pi2 / omega2) * gsl_sf_gamma((1-pi2)/omega2) * gsl_sf_gamma(n-x+(1-pi2) / omega2) * gsl_sf_gamma(n + 1/omega2);
	double ans = eta * a + (1.0 - eta) * b;
	return log(ans);
}

/*target distribution: uniform in the unit cube (side length=10)*/
double target_logdensity(void* ignore, gsl_vector* x) {
	double eta = gsl_vector_get(x, 0);
	double pi1 = gsl_vector_get(x, 1);
	double pi2 = gsl_vector_get(x, 2);
	double gamma = gsl_vector_get(x, 3);
	double ans = 0.0;
	for(int i=0; i<NR; i++)
		ans += f(gsl_vector_get(xx,i), gsl_vector_get(nn,i), eta, pi1, pi2, gamma);
	return ans;
}

int main(int argc, char** argv) {
	gsl_matrix* data = gsl_matrix_alloc(NR, 2);
	FILE* fdata = fopen(DATA_FNAME, "r");
	if(!fdata)
		return(1);
	gsl_matrix_fscanf(fdata, data);

	gsl_vector_view cv = gsl_matrix_column(data, 0);
	xx = &(cv.vector);
	cv = gsl_matrix_column(data, 1);
	nn = &(cv.vector);

	fclose(fdata);
}
