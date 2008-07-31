#ifndef __LPDF_GSL_CONT_H__
#define __LPDF_GSL_CONT_H__

/* INTERNAL UTILITY MACROS*/
#define TYPE_PAR(prefix) prefix ## _lpdf_p

#define DECLARE_2PAR(prefix, par1, par2) \
typedef struct {\
	double* par1;\
	double* par2;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(double* par1, double* par2);\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p);\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p);

#define DECLARE_1PAR(prefix, par1) \
typedef struct {\
	double* par1;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(double* par1);\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p);\
double mcmclib_ ## prefix ## _lpdf(gsl_vector* x, void* in_p);
/*END OF INTERNAL UTILITY MACROS*/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

DECLARE_1PAR(gaussian, sd)
DECLARE_1PAR(exponential, mean)
DECLARE_1PAR(laplace, a)
DECLARE_2PAR(exppow, a, b)
DECLARE_1PAR(cauchy, a)
DECLARE_1PAR(rayleigh, sigma)
DECLARE_2PAR(rayleigh_tail, a, sigma)
DECLARE_2PAR(gamma, a, b)
DECLARE_2PAR(flat, a, b)
DECLARE_2PAR(lognormal, zeta, sigma)
DECLARE_1PAR(chisq, nu)
DECLARE_2PAR(fdist, nu1, nu2)
DECLARE_1PAR(tdist, nu)
DECLARE_2PAR(beta, a, b)
DECLARE_1PAR(logistic, a)
DECLARE_2PAR(pareto, a, b)
DECLARE_2PAR(weibull, a, b)
DECLARE_2PAR(gumbel1, a, b)
DECLARE_2PAR(gumbel2, a, b)

#endif
