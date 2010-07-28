#include <assert.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "vector_queue.h"

typedef struct vector_queue_t {
  size_t dim;
  size_t max_size;
  size_t size;
  gsl_matrix* X;
  size_t *pos;
  size_t next_free;
} vector_queue;

vector_queue* vector_queue_alloc(const size_t dim, const size_t max_size) {
  vector_queue* a = (vector_queue*) malloc(sizeof(vector_queue));
  a->dim = dim;
  a->max_size = max_size;
  a->size = 0;
  a->X = gsl_matrix_alloc(max_size, dim);
  a->pos = (size_t*) malloc(max_size * sizeof(size_t));
  for(size_t i = 0; i < max_size; i++) {
    a->pos[i] = i;
  }
  a->next_free = 0;
  return a;
}

void vector_queue_free(vector_queue_t* q) {
  if(!q) return;
  gsl_matrix_free(q->X);
  free(q->pos);
  free(q);
}

int vector_queue_append(vector_queue* q, const gsl_vector* ix) {
  assert(ix->size == q->dim);
  if(q->size == q->max_size) {
    vector_queue_remove(q);
  }
  q->size = q->size + 1;
  gsl_vector_view x_v = gsl_matrix_row(q->X, q->next_free);
  gsl_vector* x = &(x_v.vector);
  gsl_vector_memcpy(x, ix);
  return 0;
}

void vector_queue_remove(vector_queue* q) {
  q->next_free = (q->next_free == 0) ? q->max_size : q->next_free - 1;
}

size_t vector_queue_size(const vector_queue* q) {
  return q->size;
}

void vector_queue_get(const vector_queue_t* q, const size_t i, gsl_vector* x) {
  assert(i < q->max_size);
  gsl_vector_view y_v = gsl_matrix_row(q->X, q->pos[i]);
  gsl_vector_memcpy(x, &(y_v.vector));
}
