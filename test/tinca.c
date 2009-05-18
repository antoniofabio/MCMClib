/**Test INCA object*/
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <inca.h>

#define v0(x) gsl_vector_get(x, 0)
#define x0 v0(x)

#define TOL 1e-6
static int check_dequal(double a, double b) {
  return (fabs(a-b) < TOL);
}

static double dtarget(void* ignore, gsl_vector* x) {
  if((x0 >= 0.0) && (x0 <= 1.0))
    return log(1.0);
  if((x0 < -10.0))
    return log(0.5);
  return log(0.0);
}

static double qd(void* ignore, gsl_vector* x, gsl_vector* y) {
  return 0.0;
}

static void sampler(mcmclib_mh_q* q, gsl_vector* x) {
  double* o = (double*) q->gamma;
  gsl_vector_set(x, 0, x0 + (*o));
}

static void update_gamma(void* in_p) {
  mcmclib_amh* p = (mcmclib_amh*) in_p;
  double *s = (double*) p->suff;
  (*s)+= v0(p->mh->x);
}

int main(int argc, char** argv) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector* x[2];
  for(int m=0; m<2; m++) {
    x[m] = gsl_vector_alloc(1);
    gsl_vector_set(x[m], 0, 0.1 * ((double) m));
  }
  double inc = 0.2;
  double suff = 0.0;
  mcmclib_mh_q* q = mcmclib_mh_q_alloc(r, sampler, &inc, qd, NULL, &inc);
  mcmclib_mh* mh = mcmclib_mh_alloc(r, dtarget, NULL, q, gsl_vector_alloc(1));
  mcmclib_amh* amh = mcmclib_amh_alloc(mh, &suff, update_gamma);
  mcmclib_inca* s = mcmclib_inca_alloc(amh, x, 2);
  double suff_check = 0.0;
  for(int n=0; n<2; n++) {
    mcmclib_inca_update(s);
    for(int m=0; m<2; m++)
      suff_check += v0(x[m]);
    assert(check_dequal(suff, suff_check));
  }
  assert(check_dequal(v0(x[0]), inc * 2));
  assert(check_dequal(v0(x[1]), inc * 2 + 0.1));

  mcmclib_inca_free(s);
  mcmclib_amh_free(amh);
  mcmclib_mh_free(mh);
  mcmclib_mh_q_free(q);
  for(int m=0; m<2; m++)
    gsl_vector_free(x[m]);
  gsl_rng_free(r);

  return 0;
}
