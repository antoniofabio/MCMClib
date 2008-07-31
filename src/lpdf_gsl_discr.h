#ifndef __LPDF_GSL_CONT_H__
#define __LPDF_GSL_CONT_H__

/* INTERNAL UTILITY MACROS*/
#define TYPE_PAR(prefix) prefix ## _lpdf_p

#define DECLARE_2PAR(prefix, type1, par1, type2, par2) \
typedef struct {\
	type1 * par1;\
	type2 * par2;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(type1 * par1, type2 * par2);\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p);\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p);

#define DECLARE_1PAR(prefix, type1, par1) \
typedef struct {\
	type1 * par1;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(type1 * par1);\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p);\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p);
/*END OF INTERNAL UTILITY MACROS*/

DECLARE_1PAR(poisson, double, mu)
DECLARE_1PAR(bernoulli, double, p)
DECLARE_2PAR(binomial, double, p, int, n)
DECLARE_2PAR(negative_binomial, double, p, int, n)
DECLARE_2PAR(pascal, double, p, int, n)
DECLARE_1PAR(geometric, double, p)
DECLARE_1PAR(logarithmic, double, p)

#endif
