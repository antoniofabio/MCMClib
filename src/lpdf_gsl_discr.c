#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include "lpdf_gsl_discr.h"

#define IMPLEMENT_1PAR_ALLOCFREE(prefix, type1, par1) \
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(type1 * par1){\
	TYPE_PAR(prefix)* ans = ( TYPE_PAR(prefix)* ) malloc(sizeof( TYPE_PAR(prefix)* ));\
	ans->par1 = par1;\
	return ans;\
}\
\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p){\
	free(p);\
}

#define IMPLEMENT_1PAR(prefix, type1, par1) \
IMPLEMENT_1PAR_ALLOCFREE(prefix, type1, par1)\
\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p) {\
	TYPE_PAR(prefix)* __p = ( TYPE_PAR(prefix)* ) in_p;\
	double ans = 0;\
	type1 par1 = *(__p->par1);\
	for(int i=0; i< x->size; i++)\
		ans += log( gsl_ran_ ## prefix ## _pdf((int) gsl_vector_get(x, i), par1) );\
	return ans;\
}

#define IMPLEMENT_2PAR_ALLOCFREE(prefix, type1, par1, type2, par2) \
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(type1 * par1, type2 * par2){\
	TYPE_PAR(prefix)* ans = ( TYPE_PAR(prefix)* ) malloc(sizeof( TYPE_PAR(prefix)* ));\
	ans->par1 = par1;\
	ans->par2 = par2;\
	return ans;\
}\
\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p){\
	free(p);\
}

#define IMPLEMENT_2PAR(prefix, type1, par1, type2, par2) \
IMPLEMENT_2PAR_ALLOCFREE(prefix, type1, par1, type2, par2)\
\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p) {\
	TYPE_PAR(prefix)* __p = ( TYPE_PAR(prefix)* ) in_p;\
	double ans = 0;\
	type1 par1 = *(__p->par1);\
	type2 par2 = *(__p->par2);\
	for(int i=0; i< x->size; i++)\
		ans += log( gsl_ran_ ## prefix ## _pdf((int) gsl_vector_get(x, i), par1, par2) );\
	return ans;\
}

IMPLEMENT_1PAR(poisson, double, mu)
IMPLEMENT_1PAR(bernoulli, double, p)
IMPLEMENT_2PAR(binomial, double, p, int, n)
IMPLEMENT_2PAR(negative_binomial, double, p, int, n)
IMPLEMENT_2PAR(pascal, double, p, int, n)
IMPLEMENT_1PAR(geometric, double, p)
IMPLEMENT_1PAR(logarithmic, double, p)
