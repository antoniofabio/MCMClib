#ifndef __MCMCLIB_COMMON_H__
#define __MCMCLIB_COMMON_H__

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/** pointer to a distribution function
*/
typedef double (*distrfun_p) (void* data, gsl_vector* x);

#endif
