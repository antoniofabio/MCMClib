/**Test hierarchical distribution*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lpdf_hierarchical.h>

#define TOL 1e-6

#define DIM 10
#define P DIM

#define x0 gsl_vector_get(x, 0)
gsl_vector *v_parent, *v_child, *x_g;
mcmclib_post_lpdf* lpost_data;

static double lprior(void* ignore, const gsl_vector* x) {
  ignore = NULL;
  if((x0 < 1.0) && (x0 > 0.0))
    return 1.0;
  else
    return 0.0;
}
static double llik(void* ignore, const gsl_vector* x) {
  ignore = NULL;
  return (x0 < gsl_vector_get(v_parent, 0)) ? 2.0 : 0.0;
}

static double lpost(double in_x) {
  gsl_vector_set(x_g, 0, in_x);
  return mcmclib_post_lpdf_compute(lpost_data, x_g);
}

int main() {
  v_parent = gsl_vector_alloc(1);
  gsl_vector_set(v_parent, 0, 0.5);
  v_child = gsl_vector_alloc(1);
  gsl_vector_set(v_child, 0, 0.4);
  int ignore[1];
  ignore[0] = 0;

  lpost_data = mcmclib_post_lpdf_alloc(v_parent, lprior, NULL,
				       llik, &v_child, (void**) &ignore, 1);

  x_g = gsl_vector_alloc(1);

  assert(lpost(0.5) == 3.0);
  assert(lpost(1.0) == 2.0);
  gsl_vector_set(v_child, 0, 0.5);
  assert(lpost(0.5) == 1.0);
  gsl_vector_set(v_child, 0, 1.5);
  assert(lpost(0.5) == 1.0);
  assert(lpost(1.5) == 0.0);
  assert(lpost(2.5) == 2.0);

  gsl_vector_free(x_g);
  mcmclib_post_lpdf_free(lpost_data);
  gsl_vector_free(v_parent);
  gsl_vector_free(v_child);
}
