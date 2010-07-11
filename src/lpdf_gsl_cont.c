/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include "lpdf_gsl_cont.h"

#define IMPLEMENT_1PAR_ALLOCFREE(prefix, par1) \
TYPE_PAR(prefix)* TYPE_METHOD(prefix, alloc)(double* par1){\
	TYPE_PAR(prefix)* ans = ( TYPE_PAR(prefix)* ) malloc(sizeof( TYPE_PAR(prefix) ));\
	ans->par1 = par1;\
	return ans;\
}\
\
void TYPE_METHOD(prefix, free)(TYPE_PAR(prefix)* p){\
	free(p);\
}

#define IMPLEMENT_1PAR(prefix, par1) \
IMPLEMENT_1PAR_ALLOCFREE(prefix, par1)\
\
double TYPE_METHOD(prefix, compute)(void* in_p, const gsl_vector* x) {\
	TYPE_PAR(prefix)* p = ( TYPE_PAR(prefix)* ) in_p;\
	double ans = 0;\
	double par1 = *(p->par1);\
	for(size_t i=0; i < x->size; i++)\
		ans += log( gsl_ran_ ## prefix ## _pdf(gsl_vector_get(x, i), par1) );\
	return ans;\
}

#define IMPLEMENT_2PAR_ALLOCFREE(prefix, par1, par2) \
TYPE_PAR(prefix)* TYPE_METHOD(prefix, alloc)(double* par1, double* par2){\
	TYPE_PAR(prefix)* ans = ( TYPE_PAR(prefix)* ) malloc(sizeof( TYPE_PAR(prefix) ));\
	ans->par1 = par1;\
	ans->par2 = par2;\
	return ans;\
}\
\
void TYPE_METHOD(prefix, free) (TYPE_PAR(prefix)* p){\
	free(p);\
}

#define IMPLEMENT_2PAR(prefix, par1, par2) \
IMPLEMENT_2PAR_ALLOCFREE(prefix, par1, par2)\
\
double TYPE_METHOD(prefix, compute)(void* in_p, const gsl_vector* x) {\
	TYPE_PAR(prefix)* p = ( TYPE_PAR(prefix)* ) in_p;\
	double ans = 0;\
	double par1 = *(p->par1);\
	double par2 = *(p->par2);\
	for(size_t i=0; i < x->size; i++)\
		ans += log( gsl_ran_ ## prefix ## _pdf(gsl_vector_get(x, i), par1, par2) );\
	return ans;\
}

IMPLEMENT_1PAR(gaussian, sd)
IMPLEMENT_1PAR(exponential, mean)
IMPLEMENT_1PAR(laplace, a)
IMPLEMENT_2PAR(exppow, a, b)
IMPLEMENT_1PAR(cauchy, a)
IMPLEMENT_1PAR(rayleigh, sigma)
IMPLEMENT_2PAR(rayleigh_tail, a, sigma)
IMPLEMENT_2PAR(gamma, a, b)
IMPLEMENT_2PAR(flat, a, b)
IMPLEMENT_2PAR(lognormal, zeta, sigma)
IMPLEMENT_1PAR(chisq, nu)
IMPLEMENT_2PAR(fdist, nu1, nu2)
IMPLEMENT_1PAR(tdist, nu)
IMPLEMENT_2PAR(beta, a, b)
IMPLEMENT_1PAR(logistic, a)
IMPLEMENT_2PAR(pareto, a, b)
IMPLEMENT_2PAR(weibull, a, b)
IMPLEMENT_2PAR(gumbel1, a, b)
IMPLEMENT_2PAR(gumbel2, a, b)
