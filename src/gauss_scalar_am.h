#ifndef __GAUSS_SCALAR_AM_H__
#define __GAUSS_SCALAR_AM_H__

#include "amh.h"

mcmclib_amh* mcmclib_gauss_scalar_am_alloc(gsl_rng* r, distrfun_p logdistr, void* logdistr_data,
					   gsl_vector* x, const double scaling, const size_t N0);

#endif
