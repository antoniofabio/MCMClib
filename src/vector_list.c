/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include "vector_list.h"

vector_list* mcmclib_vector_list_alloc() {
  vector_list* ans = (vector_list*) malloc(sizeof(vector_list));
  ans->next = NULL;
  return ans;
}

vector_list* mcmclib_vector_list_last(vector_list* i) {
  while(i->next) i = i->next;
  return i;
}

vector_list* mcmclib_vector_list_append(gsl_vector* v, vector_list* element) {
  vector_list* item = (vector_list*) malloc(sizeof(vector_list));
  item->v = v;
  item->next = NULL;
  mcmclib_vector_list_last(element)->next = item;
  return item;
}

size_t mcmclib_vector_list_length(vector_list* first) {
  size_t n = 0;
  while(first) {
    first = first->next;
    n++;
  }
  return n;
}

void mcmclib_vector_list_free(vector_list* first) {
  vector_list* save;
  while(first) {
    save = first;
    first = first->next;
    gsl_vector_free(save->v);
    free(save);
  }
}

vector_list* mcmclib_vector_list_transpose(vector_list* first) {
  size_t d = first->v->size;
  size_t n = mcmclib_vector_list_length(first);
  vector_list *last, *current;
  /*alloc result object*/
  vector_list* head = mcmclib_vector_list_alloc();
  head->v = gsl_vector_alloc(n);
  last = head;
  for(size_t i=0; i<(d-1); i++)
    last = mcmclib_vector_list_append(gsl_vector_alloc(n), last);
  /*end alloc result object*/

  /*fill out values*/
  size_t i=0;
  while(first) {
    current = head;
    size_t j=0;
    while(current) {
      gsl_vector_set(current->v, i, gsl_vector_get(first->v, j));
      current = current->next;
      j++;
    }
    first = first->next;
    i++;
  }
  /*end fill out values*/

  return head;
}

gsl_matrix* mcmclib_vector_list_asmatrix(vector_list* first) {
  size_t nr = mcmclib_vector_list_length(first);
  size_t nc = first->v->size;
  gsl_matrix* ans = gsl_matrix_alloc(nr, nc);
  for(size_t i=0; i<nr; i++) {
    for(size_t j=0; j<nc; j++)
      gsl_matrix_set(ans, i, j, gsl_vector_get(first->v, j));
    first = first->next;
  }
  return ans;
}
