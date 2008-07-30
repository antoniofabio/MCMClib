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
#include "gauss_am.h"
#include "mvnorm.h"
%}

/*vector_list methods*/
struct vector_list_str {
	gsl_vector* v;
	struct vector_list_str* next;
};
typedef struct vector_list_str vector_list;
vector_list* mcmclib_vector_list_alloc();
void mcmclib_vector_list_add(gsl_vector* v, vector_list* last);
int mcmclib_vector_list_length(vector_list* first);
void mcmclib_vector_list_free(vector_list* first);

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
	mcmclib_gauss_am_data*);
