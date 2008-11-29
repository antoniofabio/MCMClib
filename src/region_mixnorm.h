#ifndef __MCMCLIB_REGION_MIXNORM_H__
#define __MCMCLIB_REGION_MIXNORM_H__

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

/** \ingroup adaptive

Region computing function based on mixture of gaussian densities
@param x point in which compute the region
@param p pointer to an 'mcmclib_mixnorm' object
@returns to which region belongs point 'x', as an integer between 0 and K-1
*/
int mcmclib_region_mixnorm_compute(gsl_vector* x, void* p);

#endif
