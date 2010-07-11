/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __LPDF_GSL_DISCR_H__
#define __LPDF_GSL_DISCR_H__

/* INTERNAL UTILITY MACROS*/
#undef TYPE_PAR
#define TYPE_PAR(prefix) mcmclib_ ## prefix ## _lpdf
#define TYPE_METHOD(prefix, method) mcmclib_ ## prefix ## _lpdf_ ## method

#undef DECLARE_2PAR
#define DECLARE_2PAR(prefix, type1, par1, type2, par2) \
typedef struct {\
	type1 * par1;\
	type2 * par2;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* TYPE_METHOD(prefix, alloc)(type1 * par1, type2 * par2);\
void TYPE_METHOD(prefix, free)(TYPE_PAR(prefix)* p);\
double TYPE_METHOD(prefix, compute)(void* in_p, const gsl_vector* x);

#undef DECLARE_1PAR
#define DECLARE_1PAR(prefix, type1, par1) \
typedef struct {\
	type1 * par1;\
} TYPE_PAR(prefix);\
TYPE_PAR(prefix)* TYPE_METHOD(prefix, alloc)(type1 * par1);\
void TYPE_METHOD(prefix, free)(TYPE_PAR(prefix)* p);\
double TYPE_METHOD(prefix, compute)(void* in_p, const gsl_vector* x);
/*END OF INTERNAL UTILITY MACROS*/

/**\addtogroup distributions
@{*/

/**\defgroup discrete Discrete
\brief Univariate discrete probability distributions

The following probability distributions have been wrapped from the GSL:
- poisson
- bernoulli
- binomial
- negative binomial
- pascal
- geometric
- logaritmic
*/
/**@}*/

DECLARE_1PAR(poisson, double, mu)
DECLARE_1PAR(bernoulli, double, p)
DECLARE_2PAR(binomial, double, p, unsigned int, n)
DECLARE_2PAR(negative_binomial, double, p, int, n)
DECLARE_2PAR(pascal, double, p, unsigned int, n)
DECLARE_1PAR(geometric, double, p)
DECLARE_1PAR(logarithmic, double, p)

#endif
