#include <math.h>
#include "sglib.h"
#include "vector_queue.h"

#define cmp(a,b) ((b)-(a))

SGLIB_DEFINE_LIST_PROTOTYPES(vl_i, cmp, next)
SGLIB_DEFINE_LIST_FUNCTIONS(vl_i, cmp, next)

static vl_i* vl_i_alloc(gsl_vector* x) {
  vl_i* a = (vl_i*) malloc(sizeof(vl_i));
  a->x = x;
  a->next = NULL;
  return a;
}

static void vl_i_free(vl_i* p) {
  free(p);
}

typedef struct vector_queue_t {
  size_t dim;
  size_t max_size;
  size_t size;
  vl_i* head;
} vector_queue;

vector_queue* vector_queue_alloc(const size_t dim, const size_t max_size) {
  vector_queue* a = (vector_queue*) malloc(sizeof(vector_queue));
  a->dim = dim;
  a->max_size = max_size;
  a->size = 0;
  a->head = NULL;
  return a;
}

void vector_queue_free(vector_queue_t* q) {
  if(!q) return;
  SGLIB_LIST_MAP_ON_ELEMENTS(vl_i, q->head, li, next, {
      gsl_vector_free(li->x);
      vl_i_free(li);
    });
  free(q);
}

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

int vector_queue_append(vector_queue* q, const gsl_vector* ix) {
  assert(ix->size == q->dim);
  gsl_vector* x = gsl_vector_alloc(q->dim);
  gsl_vector_memcpy(x, ix);
  vl_i* li = vl_i_alloc(x);
  sglib_vl_i_add(&(q->head), li);
  if(q->size == q->max_size) {
    vector_queue_remove(q);
  }
  return 0;
}

void vector_queue_remove(vector_queue* q) {
  gsl_vector_free(q->head->x);
  q->head = q->head->next;
  free(q->head);
}

size_t vector_queue_size(const vector_queue* q) {
  return q->size;
}

void vector_queue_get(const vector_queue_t* q, const size_t lag, gsl_vector* x) {
  x = (&(q->head))[lag];
}
