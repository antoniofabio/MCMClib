#ifndef __DISTRIBS_UNIV_H__
#define __DISTRIBS_UNIV_H__

/* INTERNAL UTILITY MACROS*/
#define TYPE_PAR(prefix) prefix ## _lpdf_p

#define DECLARE_2PAR(prefix, par1, par2) \
typedef struct {\
	double* par1;\
	double* par2;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(double* par1, double* par2);\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p);

#define DECLARE_1PAR(prefix, par1) \
typedef struct {\
	double* par1;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* mcmclib_ ## prefix ## _lpdf_alloc(double* par1);\
void mcmclib_ ## prefix ## _lpdf_free(TYPE_PAR(prefix)* p);
/*END OF INTERNAL UTILITY MACROS*/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

DECLARE_2PAR(gaussian, mean, precision)
DECLARE_1PAR(exponential, mean)

#endif
