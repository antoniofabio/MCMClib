%{
#include <gsl/gsl_matrix.h>
%}

/*
GSL_MATRIX
*/
typedef struct {
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
} gsl_matrix;

%array_functions(gsl_matrix*, matrixArray);

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

int gsl_matrix_fread (FILE * istream, gsl_matrix * m) ;
int gsl_matrix_fwrite (FILE * ostream, const gsl_matrix * m) ;
int gsl_matrix_fscanf (FILE * istream, gsl_matrix * m);
int gsl_matrix_fprintf (FILE * ostream, const gsl_matrix * m, const char * format);
 
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
