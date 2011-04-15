#include <assert.h>
#include "gauss_scalar_am.h"

#define VECTOR_MAP(x, op) for(size_t i = 0; i < x->size; i++) { \
  double xi = gsl_vector_get(x, i);			\
  op;							\
  gsl_vector_set(x, i, xi);				\
  }

typedef struct {
  gsl_vector* sd;
  double scaling;
} gamma_t;

static gamma_t* gamma_t_alloc(const size_t dim, const double s0, const double scaling) {
  gamma_t* g = (gamma_t*) malloc(sizeof(gamma_t));
  g->sd = gsl_vector_alloc(dim);
  gsl_vector_set_all(g->sd, s0);
  g->scaling = scaling;
  return g;
}
static void gamma_t_free(void* in_gamma) {
  if(!in_gamma) return;
  gamma_t* g = (gamma_t*) in_gamma;
  gsl_vector_free(g->sd);
  free(g);
}
static void sample(mcmclib_mh_q* q, gsl_vector* x) {
  gamma_t* g = (gamma_t*) q->gamma;
  VECTOR_MAP(x, xi += gsl_ran_gaussian(q->r, gsl_vector_get(g->sd, i) * g->scaling));
}
static mcmclib_mh_q* mh_q_alloc(gsl_rng* r, const size_t dim, const double s0, const double scaling) {
  return mcmclib_mh_q_alloc(r, sample, NULL, gamma_t_alloc(dim, s0, scaling), gamma_t_free);
}

typedef struct {
  size_t n;
  gsl_vector* sum_x;
  gsl_vector* sum_xsq;
  gsl_vector* tmp;
  gsl_vector* mean;
  gsl_vector* var;
  double scaling_factor;
} suff_t;

static void gsl_vector_square(gsl_vector* x) {
  VECTOR_MAP(x, xi *= xi);
}

static suff_t* suff_t_alloc(const size_t size) {
  suff_t* ans = (suff_t*) malloc(sizeof(suff_t));
  ans->n = 0;
  ans->sum_x = gsl_vector_alloc(size);
  gsl_vector_set_zero(ans->sum_x);
  ans->sum_xsq = gsl_vector_alloc(size);
  gsl_vector_set_zero(ans->sum_xsq);
  ans->tmp = gsl_vector_alloc(size);
  ans->mean = gsl_vector_alloc(size);
  ans->var = gsl_vector_alloc(size);
  return ans;
}

static void suff_t_free(void* in_s) {
  suff_t* s = (suff_t*) in_s;
  if(!s) return;
  gsl_vector_free(s->var);
  gsl_vector_free(s->mean);
  gsl_vector_free(s->tmp);
  gsl_vector_free(s->sum_x);
  gsl_vector_free(s->sum_xsq);
  free(s);
}

static void suff_update(void* s, const gsl_vector* x) {
  suff_t* suff = (suff_t*) s;
  suff->n++;
  gsl_vector_add(suff->sum_x, x);
  gsl_vector_memcpy(suff->tmp, x);
  gsl_vector_square(suff->tmp);
  gsl_vector_add(suff->sum_xsq, suff->tmp);
  gsl_vector_memcpy(suff->mean, suff->sum_x);
  gsl_vector_scale(suff->mean, 1.0 / (double) (suff->n));
  gsl_vector_square(suff->mean);
  gsl_vector_memcpy(suff->var, suff->sum_xsq);
  gsl_vector_scale(suff->var, 1.0 / (double) (suff->n));
  gsl_vector_sub(suff->var, suff->mean);
  assert((!gsl_vector_isneg(suff->var)));
}

static void gamma_update(void* s, void* in_g) {
  suff_t* suff = (suff_t*) s;
  gamma_t* g = (gamma_t*) in_g;
  gsl_vector_memcpy(g->sd, suff->var);
  VECTOR_MAP(g->sd, xi = fmax(1e-3, sqrt(xi)));
}

mcmclib_amh* mcmclib_gauss_scalar_am_alloc(gsl_rng* r, distrfun_p logdistr, void* logdistr_data,
					   gsl_vector* x, const double s0, const double scaling, const size_t N0) {
  mcmclib_amh* ans = mcmclib_amh_alloc(mcmclib_mh_alloc(r, logdistr, logdistr_data, mh_q_alloc(r, x->size, s0, scaling), x), N0,
				       suff_t_alloc(x->size), suff_t_free, suff_update, gamma_update);
  return ans;
}
