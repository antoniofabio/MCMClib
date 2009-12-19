%{
#include <gsl/gsl_vector.h>
%}

/*
GSL_VECTOR
*/
typedef struct {
  size_t size;
  size_t stride;
  double *data;
} gsl_vector;

%array_functions(gsl_vector*, vectorArray);

double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);

double *gsl_vector_ptr (gsl_vector * v, const size_t i);
const double *gsl_vector_const_ptr (const gsl_vector * v, const size_t i);

void gsl_vector_set_zero (gsl_vector * v);
void gsl_vector_set_all (gsl_vector * v, double x);
int gsl_vector_set_basis (gsl_vector * v, size_t i);

int gsl_vector_fread (FILE * istream, gsl_vector * v);
int gsl_vector_fwrite (FILE * ostream, const gsl_vector * v);
int gsl_vector_fscanf (FILE * istream, gsl_vector * v);
int gsl_vector_fprintf (FILE * ostream, const gsl_vector * v,
			const char *format);

int gsl_vector_memcpy (gsl_vector * dest, const gsl_vector * src);

int gsl_vector_reverse (gsl_vector * v);

int gsl_vector_swap (gsl_vector * v, gsl_vector * w);
int gsl_vector_swap_elements (gsl_vector * v, const size_t i, const size_t j);

double gsl_vector_max (const gsl_vector * v);
double gsl_vector_min (const gsl_vector * v);
void gsl_vector_minmax (const gsl_vector * v, double * min_out, double * max_out);

size_t gsl_vector_max_index (const gsl_vector * v);
size_t gsl_vector_min_index (const gsl_vector * v);
void gsl_vector_minmax_index (const gsl_vector * v, size_t * imin, size_t * imax);

int gsl_vector_add (gsl_vector * a, const gsl_vector * b);
int gsl_vector_sub (gsl_vector * a, const gsl_vector * b);
int gsl_vector_mul (gsl_vector * a, const gsl_vector * b);
int gsl_vector_div (gsl_vector * a, const gsl_vector * b);
int gsl_vector_scale (gsl_vector * a, const double x);
int gsl_vector_add_constant (gsl_vector * a, const double x);

int gsl_vector_isnull (const gsl_vector * v);
int gsl_vector_ispos (const gsl_vector * v);
int gsl_vector_isneg (const gsl_vector * v);
int gsl_vector_isnonneg (const gsl_vector * v);

%extend gsl_vector {
  gsl_vector(size_t n) {
    return gsl_vector_alloc(n);
  }
  ~gsl_vector() {
    gsl_vector_free($self);
  }
};
