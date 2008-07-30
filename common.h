#ifndef __MCMCLIB_COMMON_H__
#define __MCMCLIB_COMMON_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>

#define isfinite(x) !(isnan((x)) || isinf((x)))

#endif
