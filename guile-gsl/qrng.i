%{
#include <gsl/gsl_qrng.h>
%}

/*QRNG: Quasi-random sequences*/
typedef struct {} gsl_qrng;
%extend gsl_qrng {
  gsl_qrng(gsl_qrng_type* T, unsigned int d) {
    return gsl_qrng_alloc(T, d);
  }
  ~gsl_qrng() {
    gsl_qrng_free($self);
  }
  int get_value(gsl_vector* x) {
    return gsl_qrng_get($self, x->data);
  }
}
void gsl_qrng_init (gsl_qrng * q);
const char * gsl_qrng_name (const gsl_qrng * q);
size_t gsl_qrng_size (const gsl_qrng * q);
void * gsl_qrng_state (const gsl_qrng * q);
int gsl_qrng_memcpy (gsl_qrng * dest, const gsl_qrng * src);
%newobject gsl_qrng_clone;
gsl_qrng * gsl_qrng_clone (const gsl_qrng * q);
const gsl_qrng_type* gsl_qrng_niederreiter_2;
const gsl_qrng_type* gsl_qrng_sobol;
const gsl_qrng_type* gsl_qrng_halton;
const gsl_qrng_type* gsl_qrng_reversehalton;
