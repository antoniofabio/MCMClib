%module gsl

%scheme %{(load-extension "libguilegsl.so" "scm_init_swig_gsl_module")%}

%{
#include <string.h>
#include <libguile.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

static void guile_gsl_err_handler(const char * reason,
				  const char * file,
				  int line,
				  int gsl_errno) {
  static char msg[2048];
  sprintf(msg, "%s:%d", file, line);
  scm_error_scm(scm_misc_error_key,
		scm_from_locale_string(msg),
		scm_from_locale_string(reason), SCM_BOOL_F, SCM_BOOL_F);
}
%}

%init %{
  gsl_set_error_handler(guile_gsl_err_handler);
%}

FILE* fopen(const char*, const char*);
void fclose(FILE*);

/*
GSL_VECTOR
*/
typedef struct {
  size_t size;
  size_t stride;
  double *data;
} gsl_vector;

double gsl_vector_get (const gsl_vector * v, const size_t i);
void gsl_vector_set (gsl_vector * v, const size_t i, double x);

double *gsl_vector_ptr (gsl_vector * v, const size_t i);
const double *gsl_vector_const_ptr (const gsl_vector * v, const size_t i);

void gsl_vector_set_zero (gsl_vector * v);
void gsl_vector_set_all (gsl_vector * v, double x);
int gsl_vector_set_basis (gsl_vector * v, size_t i);

int gsl_vector_fread (FILE * stream, gsl_vector * v);
int gsl_vector_fwrite (FILE * stream, const gsl_vector * v);
int gsl_vector_fscanf (FILE * stream, gsl_vector * v);
int gsl_vector_fprintf (FILE * stream, const gsl_vector * v,
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

/*
GSL_MATRIX
*/
typedef struct {
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
} gsl_matrix;

gsl_matrix * 
gsl_matrix_calloc (const size_t n1, const size_t n2);

gsl_matrix * 
gsl_matrix_alloc_from_block (gsl_block * b, 
                                   const size_t offset, 
                                   const size_t n1, 
                                   const size_t n2, 
                                   const size_t d2);

gsl_matrix * 
gsl_matrix_alloc_from_matrix (gsl_matrix * m,
                                    const size_t k1, 
                                    const size_t k2,
                                    const size_t n1, 
                                    const size_t n2);

gsl_vector * 
gsl_vector_alloc_row_from_matrix (gsl_matrix * m,
                                        const size_t i);

gsl_vector * 
gsl_vector_alloc_col_from_matrix (gsl_matrix * m,
                                        const size_t j);

/* Views */

_gsl_matrix_view 
gsl_matrix_submatrix (gsl_matrix * m, 
                            const size_t i, const size_t j, 
                            const size_t n1, const size_t n2);

_gsl_vector_view 
gsl_matrix_row (gsl_matrix * m, const size_t i);

_gsl_vector_view 
gsl_matrix_column (gsl_matrix * m, const size_t j);

_gsl_vector_view 
gsl_matrix_diagonal (gsl_matrix * m);

_gsl_vector_view 
gsl_matrix_subdiagonal (gsl_matrix * m, const size_t k);

_gsl_vector_view 
gsl_matrix_superdiagonal (gsl_matrix * m, const size_t k);

_gsl_vector_view
gsl_matrix_subrow (gsl_matrix * m, const size_t i,
                         const size_t offset, const size_t n);

_gsl_vector_view
gsl_matrix_subcolumn (gsl_matrix * m, const size_t j,
                            const size_t offset, const size_t n);

_gsl_matrix_view
gsl_matrix_view_array (double * base,
                             const size_t n1, 
                             const size_t n2);

_gsl_matrix_view
gsl_matrix_view_array_with_tda (double * base, 
                                      const size_t n1, 
                                      const size_t n2,
                                      const size_t tda);


_gsl_matrix_view
gsl_matrix_view_vector (gsl_vector * v,
                              const size_t n1, 
                              const size_t n2);

_gsl_matrix_view
gsl_matrix_view_vector_with_tda (gsl_vector * v,
                                       const size_t n1, 
                                       const size_t n2,
                                       const size_t tda);


_gsl_matrix_const_view 
gsl_matrix_const_submatrix (const gsl_matrix * m, 
                                  const size_t i, const size_t j, 
                                  const size_t n1, const size_t n2);

_gsl_vector_const_view 
gsl_matrix_const_row (const gsl_matrix * m, 
                            const size_t i);

_gsl_vector_const_view 
gsl_matrix_const_column (const gsl_matrix * m, 
                               const size_t j);

_gsl_vector_const_view
gsl_matrix_const_diagonal (const gsl_matrix * m);

_gsl_vector_const_view 
gsl_matrix_const_subdiagonal (const gsl_matrix * m, 
                                    const size_t k);

_gsl_vector_const_view 
gsl_matrix_const_superdiagonal (const gsl_matrix * m, 
                                      const size_t k);

_gsl_vector_const_view
gsl_matrix_const_subrow (const gsl_matrix * m, const size_t i,
                               const size_t offset, const size_t n);

_gsl_vector_const_view
gsl_matrix_const_subcolumn (const gsl_matrix * m, const size_t j,
                                  const size_t offset, const size_t n);

_gsl_matrix_const_view
gsl_matrix_const_view_array (const double * base,
                                   const size_t n1, 
                                   const size_t n2);

_gsl_matrix_const_view
gsl_matrix_const_view_array_with_tda (const double * base, 
                                            const size_t n1, 
                                            const size_t n2,
                                            const size_t tda);

_gsl_matrix_const_view
gsl_matrix_const_view_vector (const gsl_vector * v,
                                    const size_t n1, 
                                    const size_t n2);

_gsl_matrix_const_view
gsl_matrix_const_view_vector_with_tda (const gsl_vector * v,
                                             const size_t n1, 
                                             const size_t n2,
                                             const size_t tda);

/* Operations */

double   gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j);
void    gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x);

double * gsl_matrix_ptr(gsl_matrix * m, const size_t i, const size_t j);
const double * gsl_matrix_const_ptr(const gsl_matrix * m, const size_t i, const size_t j);

void gsl_matrix_set_zero (gsl_matrix * m);
void gsl_matrix_set_identity (gsl_matrix * m);
void gsl_matrix_set_all (gsl_matrix * m, double x);

int gsl_matrix_fread (FILE * stream, gsl_matrix * m) ;
int gsl_matrix_fwrite (FILE * stream, const gsl_matrix * m) ;
int gsl_matrix_fscanf (FILE * stream, gsl_matrix * m);
int gsl_matrix_fprintf (FILE * stream, const gsl_matrix * m, const char * format);
 
int gsl_matrix_memcpy(gsl_matrix * dest, const gsl_matrix * src);
int gsl_matrix_swap(gsl_matrix * m1, gsl_matrix * m2);

int gsl_matrix_swap_rows(gsl_matrix * m, const size_t i, const size_t j);
int gsl_matrix_swap_columns(gsl_matrix * m, const size_t i, const size_t j);
int gsl_matrix_swap_rowcol(gsl_matrix * m, const size_t i, const size_t j);
int gsl_matrix_transpose (gsl_matrix * m);
int gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src);

double gsl_matrix_max (const gsl_matrix * m);
double gsl_matrix_min (const gsl_matrix * m);
void gsl_matrix_minmax (const gsl_matrix * m, double * min_out, double * max_out);

void gsl_matrix_max_index (const gsl_matrix * m, size_t * imax, size_t *jmax);
void gsl_matrix_min_index (const gsl_matrix * m, size_t * imin, size_t *jmin);
void gsl_matrix_minmax_index (const gsl_matrix * m, size_t * imin, size_t * jmin, size_t * imax, size_t * jmax);

int gsl_matrix_isnull (const gsl_matrix * m);
int gsl_matrix_ispos (const gsl_matrix * m);
int gsl_matrix_isneg (const gsl_matrix * m);
int gsl_matrix_isnonneg (const gsl_matrix * m);

int gsl_matrix_add (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_sub (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_mul_elements (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_div_elements (gsl_matrix * a, const gsl_matrix * b);
int gsl_matrix_scale (gsl_matrix * a, const double x);
int gsl_matrix_add_constant (gsl_matrix * a, const double x);
int gsl_matrix_add_diagonal (gsl_matrix * a, const double x);

%extend gsl_matrix {
  gsl_matrix(size_t n1, size_t n2) {
    return gsl_matrix_alloc(n1, n2);
  }
  ~gsl_matrix() {
    gsl_matrix_free($self);
  }
};

/*
RNG
*/
const gsl_rng_type *gsl_rng_borosh13;
const gsl_rng_type *gsl_rng_coveyou;
const gsl_rng_type *gsl_rng_cmrg;
const gsl_rng_type *gsl_rng_fishman18;
const gsl_rng_type *gsl_rng_fishman20;
const gsl_rng_type *gsl_rng_fishman2x;
const gsl_rng_type *gsl_rng_gfsr4;
const gsl_rng_type *gsl_rng_knuthran;
const gsl_rng_type *gsl_rng_knuthran2;
const gsl_rng_type *gsl_rng_knuthran2002;
const gsl_rng_type *gsl_rng_lecuyer21;
const gsl_rng_type *gsl_rng_minstd;
const gsl_rng_type *gsl_rng_mrg;
const gsl_rng_type *gsl_rng_mt19937;
const gsl_rng_type *gsl_rng_mt19937_1999;
const gsl_rng_type *gsl_rng_mt19937_1998;
const gsl_rng_type *gsl_rng_r250;
const gsl_rng_type *gsl_rng_ran0;
const gsl_rng_type *gsl_rng_ran1;
const gsl_rng_type *gsl_rng_ran2;
const gsl_rng_type *gsl_rng_ran3;
const gsl_rng_type *gsl_rng_rand;
const gsl_rng_type *gsl_rng_rand48;
const gsl_rng_type *gsl_rng_random128_bsd;
const gsl_rng_type *gsl_rng_random128_glibc2;
const gsl_rng_type *gsl_rng_random128_libc5;
const gsl_rng_type *gsl_rng_random256_bsd;
const gsl_rng_type *gsl_rng_random256_glibc2;
const gsl_rng_type *gsl_rng_random256_libc5;
const gsl_rng_type *gsl_rng_random32_bsd;
const gsl_rng_type *gsl_rng_random32_glibc2;
const gsl_rng_type *gsl_rng_random32_libc5;
const gsl_rng_type *gsl_rng_random64_bsd;
const gsl_rng_type *gsl_rng_random64_glibc2;
const gsl_rng_type *gsl_rng_random64_libc5;
const gsl_rng_type *gsl_rng_random8_bsd;
const gsl_rng_type *gsl_rng_random8_glibc2;
const gsl_rng_type *gsl_rng_random8_libc5;
const gsl_rng_type *gsl_rng_random_bsd;
const gsl_rng_type *gsl_rng_random_glibc2;
const gsl_rng_type *gsl_rng_random_libc5;
const gsl_rng_type *gsl_rng_randu;
const gsl_rng_type *gsl_rng_ranf;
const gsl_rng_type *gsl_rng_ranlux;
const gsl_rng_type *gsl_rng_ranlux389;
const gsl_rng_type *gsl_rng_ranlxd1;
const gsl_rng_type *gsl_rng_ranlxd2;
const gsl_rng_type *gsl_rng_ranlxs0;
const gsl_rng_type *gsl_rng_ranlxs1;
const gsl_rng_type *gsl_rng_ranlxs2;
const gsl_rng_type *gsl_rng_ranmar;
const gsl_rng_type *gsl_rng_slatec;
const gsl_rng_type *gsl_rng_taus;
const gsl_rng_type *gsl_rng_taus2;
const gsl_rng_type *gsl_rng_taus113;
const gsl_rng_type *gsl_rng_transputer;
const gsl_rng_type *gsl_rng_tt800;
const gsl_rng_type *gsl_rng_uni;
const gsl_rng_type *gsl_rng_uni32;
const gsl_rng_type *gsl_rng_vax;
const gsl_rng_type *gsl_rng_waterman14;
const gsl_rng_type *gsl_rng_zuf;

const gsl_rng_type ** gsl_rng_types_setup(void);

const gsl_rng_type *gsl_rng_default;

typedef struct {} gsl_rng;

%extend gsl_rng {
  gsl_rng(gsl_rng_type* type) {
    return gsl_rng_alloc(type);
  }
  ~gsl_rng() {
    gsl_rng_free($self);
  }
}
