#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include "distribs_univ.h"

#define IMPLEMENT_1PAR_ALLOCFREE(prefix, par1) \
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(double* par1){\
	TYPE_PAR(prefix)* ans = ( TYPE_PAR(prefix)* ) malloc(sizeof( TYPE_PAR(prefix)* ));\
	ans->par1 = par1;\
	return ans;\
}\
\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p){\
	free(p);\
}

#define IMPLEMENT_1PAR(prefix, par1) \
IMPLEMENT_1PAR_ALLOCFREE(prefix, par1)\
\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p) {\
	TYPE_PAR(prefix)* p = ( TYPE_PAR(prefix)* ) in_p;\
	double ans = 0;\
	double par1 = *(p->par1);\
	for(int i=0; i< x->size; i++)\
		ans += log( gsl_ran_ ## prefix ## _pdf(gsl_vector_get(x, i), par1) );\
	return ans;\
}

#define IMPLEMENT_2PAR_ALLOCFREE(prefix, par1, par2) \
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(double* par1, double* par2){\
	TYPE_PAR(prefix)* ans = ( TYPE_PAR(prefix)* ) malloc(sizeof( TYPE_PAR(prefix)* ));\
	ans->par1 = par1;\
	ans->par2 = par2;\
	return ans;\
}\
\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p){\
	free(p);\
}

/*
GAUSSIAN DISTRIBUTION
*/
IMPLEMENT_2PAR_ALLOCFREE(gaussian, mean, precision)

double mcmclib_gaussian_lpdf(gsl_vector* x, void* in_p) {
	TYPE_PAR(gaussian)* p = (TYPE_PAR(gaussian)*) in_p;
	double ans = 0;
	double mean = *(p->mean);
	double sd = sqrt(1.0 / (*(p->precision)) );
	for(int i=0; i<x->size; i++)
		ans += log( gsl_ran_gaussian_pdf(gsl_vector_get(x, i) - mean, sd) );
	return ans;
}

/*
EXPONENTIAL DISTRIBUTION
*/
IMPLEMENT_1PAR(exponential, mean)
