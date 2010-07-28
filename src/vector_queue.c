#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include "vector_queue.h"

typedef struct mcmclib_vector_queue_t {
  size_t dim;
  size_t max_size;
  size_t size;
  gsl_matrix* X;
  gsl_permutation* perm;
  size_t next_free;
} vector_queue;

static void perm_next(gsl_permutation* p) {
  const size_t size = gsl_permutation_size(p);
  size_t last = gsl_permutation_get(p, size - 1);
  for(size_t i = (size-1); i > 0; i--) {
    p->data[i] = p->data[i-1];
  }
  p->data[0] = last;
}

#if 0
static void vector_queue_fprintf(FILE* f, mcmclib_vector_queue* q) {
  fprintf(f, "permutation: ");
  gsl_permutation_fprintf(f, q->perm, "%zu ");
  fprintf(f, "\nnext_free: %zd\n", q->next_free);
  fprintf(f, "X:\n");
  gsl_matrix_fprintf(stdout, q->X, "%.3f");
}
#endif

mcmclib_vector_queue* mcmclib_vector_queue_alloc(const size_t dim, const size_t max_size) {
  vector_queue* a = (vector_queue*) malloc(sizeof(vector_queue));
  a->dim = dim;
  a->max_size = max_size;
  a->size = 0;
  a->X = gsl_matrix_alloc(max_size, dim);
  gsl_matrix_set_all(a->X, GSL_NAN);
  a->perm = gsl_permutation_alloc(max_size);
  gsl_permutation_init(a->perm);
  a->next_free = max_size - 1;
  return a;
}

void mcmclib_vector_queue_free(vector_queue* q) {
  if(!q) return;
  gsl_matrix_free(q->X);
  gsl_permutation_free(q->perm);
  free(q);
}

int mcmclib_vector_queue_append(mcmclib_vector_queue* q, const gsl_vector* ix) {
  assert(ix->size == q->dim);
  gsl_vector_view x_v = gsl_matrix_row(q->X, q->next_free);
  gsl_vector_memcpy(&(x_v.vector), ix);
  if(q->size < q->max_size) {
    q->size = q->size + 1;
  } else {
    perm_next(q->perm);
  }
  q->next_free = (q->next_free == 0) ? (q->max_size - 1) : q->next_free - 1;
  return 0;
}

size_t mcmclib_vector_queue_size(const mcmclib_vector_queue* q) {
  return q->size;
}

int mcmclib_vector_queue_get(const mcmclib_vector_queue* q, const size_t i, gsl_vector* x) {
  assert(i < q->max_size);
  if(i == q->size) {
    static char msg[524];
    sprintf(msg, "requested element %zd, but current queue size is %zd", i, q->size);
    GSL_ERROR(msg, GSL_EDOM);
  }
  gsl_vector_view y_v;
  if(q->size == q->max_size) {
    y_v = gsl_matrix_row(q->X, gsl_permutation_get(q->perm, i));
  } else {
    y_v = gsl_matrix_row(q->X, q->max_size - q->size + i);
  }
  gsl_vector_memcpy(x, &(y_v.vector));
  return 0;
}
