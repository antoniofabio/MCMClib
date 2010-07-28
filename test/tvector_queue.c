#include <vector_queue.h>
#include "CuTest.h"

#define DIM 1
#define MAX_SIZE 2

#define x0 (x->data[0])

static gsl_vector* x;
static gsl_vector* y;
static mcmclib_vector_queue* q;

static void append_value(double xi) {
  gsl_vector_set(x, 0, xi);
  mcmclib_vector_queue_append(q, x);
}

static double get_value(size_t l) {
  mcmclib_vector_queue_get(q, l, y);
  return gsl_vector_get(y, 0);
}

void Testvector_queue(CuTest* tc) {
  x = gsl_vector_alloc(1);
  y = gsl_vector_alloc(1);
  q = mcmclib_vector_queue_alloc(DIM, MAX_SIZE);
  CuAssertIntEquals(tc, 0, mcmclib_vector_queue_size(q));
  append_value(1.0);
  CuAssertIntEquals(tc, 1, mcmclib_vector_queue_size(q));
  CuAssertTrue(tc, get_value(0) == 1.0);
  append_value(2.0);
  CuAssertTrue(tc, get_value(0) == 2.0);
  CuAssertTrue(tc, get_value(1) == 1.0);
  append_value(3.0);
  CuAssertTrue(tc, get_value(0) == 3.0);
  CuAssertTrue(tc, get_value(1) == 2.0);

  mcmclib_vector_queue_free(q);
  gsl_vector_free(x);
  gsl_vector_free(y);
}
