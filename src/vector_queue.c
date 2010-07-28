#include <assert.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include "vector_queue.h"

typedef struct vector_queue_t {
  size_t dim;
  size_t max_size;
  size_t size;
  gsl_matrix* X;
  gsl_permutation* perm;
  size_t next_free;
} vector_queue;

vector_queue* vector_queue_alloc(const size_t dim, const size_t max_size) {
  vector_queue* a = (vector_queue*) malloc(sizeof(vector_queue));
  a->dim = dim;
  a->max_size = max_size;
  a->size = 0;
  a->X = gsl_matrix_alloc(max_size, dim);
  a->perm = gsl_permutation_alloc(max_size);
  gsl_permutation_init(a->perm);
  gsl_permutation_reverse(a->perm);
  a->next_free = max_size - 1;
  return a;
}

void vector_queue_free(vector_queue_t* q) {
  if(!q) return;
  gsl_matrix_free(q->X);
  gsl_permutation_free(q->perm);
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
  q->next_free = (q->next_free == 0) ? q->max_size : q->next_free - 1;
  return 0;
}

void vector_queue_remove(vector_queue* q) {
  assert(q->size > 0);
  q->size -= 1;
  q->next_free = (q->next_free + 1) % q->max_size;
}

size_t vector_queue_size(const vector_queue* q) {
  return q->size;
}

int vector_queue_get(const vector_queue_t* q, const size_t i, gsl_vector* x) {
  assert(i < q->max_size);
  if(i == q->size) {
    static char msg[524];
    sprintf(msg, "requested element %zd, but current queue size is %zd", i, q->size);
    GSL_ERROR(msg, GSL_EDOM);
  }
  gsl_vector_view y_v = gsl_matrix_row(q->X, gsl_permutation_get(q->perm, i));
  gsl_vector_memcpy(x, &(y_v.vector));
  return 0;
}
