/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2010 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __VECTOR_QUEUE_H__
#define __VECTOR_QUEUE_H__

/**\addtogroup utils
@{
\defgroup utils misc utilities
@{*/

#include <gsl/gsl_vector.h>

typedef struct vl_i {
  gsl_vector* x;
  struct vl_i* next;
} vl_i;

typedef struct vector_queue_t vector_queue_t;

vector_queue_t* vector_queue_alloc(const size_t dim, const size_t max_size);
void vector_queue_free(vector_queue_t* q);
int vector_queue_append(vector_queue_t* q, const gsl_vector* x);
void vector_queue_remove(vector_queue_t* q);
size_t vector_queue_size(const vector_queue_t* q);
#define VECTOR_QUEUE_MAP(q, elt, op) SGLIB_LIST_MAP_ON_ELEMENTS(vl_i, q->head, elt, next, op)

/**@}*/
/**@}*/

#endif
