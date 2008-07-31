%module mcmclib

%include "typemaps.i"
%include "carrays.i"

%array_functions(double, doubleArray);
%array_functions(char *, stringArray);
%array_functions(int, intArray);

%init %{
	if (Tcl_InitStubs(interp, "8.4", 0) == NULL) {
		return TCL_ERROR;
	}
%}

%{
/* Put header files here or function declarations like below */
#include "common.h"
#include "vector_list.h"
#include "vector_stats.h"
#include "gauss_rw.h"
#include "gauss_am.h"
#include "mvnorm.h"
#include "lpdf_gsl_cont.h"
#include "lpdf_gsl_discr.h"

double mcmclib_test_dunif(gsl_vector* px, const void* ignore) {
	for(int i=0; i<px->size; i++) {
		double x = gsl_vector_get(px, i);
		if((x < 0.0) || (x > 1.0))
			return log(0);
	}
	return 0;
}

%}

%constant double mcmclib_test_dunif(gsl_vector*, const void*);

/*vector_list methods*/
struct vector_list_str {
	gsl_vector* v;
	struct vector_list_str* next;
};
typedef struct vector_list_str vector_list;
vector_list* mcmclib_vector_list_alloc();
vector_list* mcmclib_vector_list_last(vector_list* i);
vector_list* mcmclib_vector_list_append(gsl_vector* v, vector_list* last);
int mcmclib_vector_list_length(vector_list* first);
void mcmclib_vector_list_free(vector_list* first);
vector_list* mcmclib_vector_list_transpose(vector_list* first);
gsl_matrix* mcmclib_vector_list_asmatrix(vector_list* first);

/*vectorial statistics functions*/
void mcmclib_matrix_colmeans(gsl_matrix* m, gsl_vector* out);
void mcmclib_matrix_rowmeans(gsl_matrix* m, gsl_vector* out);
void mcmclib_matrix_covariance(gsl_matrix* m, gsl_matrix* out);
void mcmclib_covariance_update(gsl_matrix* cov, gsl_vector* mean, int* n, gsl_vector* x);

/*multivariate gaussian variates*/
void mcmclib_mvnorm(const gsl_rng* r, const gsl_matrix* sigma, gsl_vector* out);
void mcmclib_mvnorm_chol(const gsl_rng* r, const gsl_matrix* sigma_chol, gsl_vector* out);

/*Gaussian random walk*/
mcmclib_gauss_rw_data* mcmclib_gauss_rw_alloc(double step_size, int dim);
void mcmclib_gauss_rw_free(mcmclib_gauss_rw_data* p);
int mcmclib_gauss_rw(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	mcmclib_gauss_rw_data* e);

/*Adaptive Metropolis algorithm*/
typedef struct {
	gsl_matrix* sigma_zero;
	int t0;
	gsl_vector* mean;
	gsl_matrix* cov;
	int t;
	gsl_vector* old;
} mcmclib_gauss_am_data;
mcmclib_gauss_am_data* mcmclib_gauss_am_alloc(const gsl_matrix* sigma_zero, int t0);
void mcmclib_gauss_am_free(mcmclib_gauss_am_data* p);
int mcmclib_gauss_am(const gsl_rng* r,
	double (*loglik) (gsl_vector* x, const void* data), gsl_vector* x, const void* data,
	mcmclib_gauss_am_data* e);

%include "gsl_lpdf_cont.i"
%include "gsl_lpdf_discr.i"
%include "gsl.i"
