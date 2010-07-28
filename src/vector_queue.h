/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_VECTOR_QUEUE_H__
#define __MCMCLIB_VECTOR_QUEUE_H__

/**\addtogroup utils
@{
\defgroup utils misc utilities
@{*/

#include <gsl/gsl_vector.h>

typedef struct mcmclib_vector_queue_t mcmclib_vector_queue;

mcmclib_vector_queue* mcmclib_vector_queue_alloc(const size_t dim, const size_t max_size);
void mcmclib_vector_queue_free(mcmclib_vector_queue* q);
int mcmclib_vector_queue_append(mcmclib_vector_queue* q, const gsl_vector* x);
size_t mcmclib_vector_queue_size(const mcmclib_vector_queue* q);
int mcmclib_vector_queue_get(const mcmclib_vector_queue* q, const size_t lag, gsl_vector* x);

/**@}*/
/**@}*/

#endif
