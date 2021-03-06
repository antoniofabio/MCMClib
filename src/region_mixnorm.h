/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009,2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_REGION_MIXNORM_H__
#define __MCMCLIB_REGION_MIXNORM_H__

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

/** \ingroup adaptive

\brief Region computing function based on mixture of gaussian densities

Intended to be used with \ref RAPT (and derived) samplers
@param x point in which compute the region
@param p pointer to an 'mcmclib_mixnorm' object
@returns to which region belongs point 'x', as an integer between 0 and K-1
*/
size_t mcmclib_region_mixnorm_compute(const gsl_vector* x, void* p);

#endif
