#include <math.h>
#include "sglib.h"
#include "vector_queue.h"

typedef struct vector_queue_t {
  size_t dim;
  size_t max_size;
  size_t size;
  gsl_vector* head;
} vector_queue;

vector_queue* vector_queue_alloc(const size_t dim, const size_t max_size) {
  vector_queue* a = (vector_queue*) malloc(sizeof(vector_queue));
  a->dim = dim;
  a->max_size = max_size;
  a->size = 0;
  a->head = NULL;
  return a;
}

int vector_queue_append(vector_queue* q, const gsl_vector* x) {
  q->size = min(q->max_size, q->size+1);
  q->head = NULL;
  return 0;
}

void vector_queue_remove(vector_queue* q) {
  q->head = NULL;
}

size_t vector_queue_size(const vector_queue* q) {
  return q->size;
}

void vector_queue_get(const vector_queue_t* q, const size_t lag, gsl_vector* x) {
  x = (&(q->head))[lag];
}
