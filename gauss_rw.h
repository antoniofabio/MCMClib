#ifndef __GAUSS_RW_H__
#define __GAUSS_RW_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#define isfinite(x) !(isnan((x)) || isinf((x)))

int mcmclib_gauss_rw(gsl_rng* r,
	double (*loglik) (gsl_vector* x, void* data), gsl_vector* x, void* data,
	double step_size);

#endif
