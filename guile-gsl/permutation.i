%{
#include <gsl/gsl_permutation.h>
%}

struct gsl_permutation_struct
{
  size_t size;
  size_t *data;
};

typedef struct gsl_permutation_struct gsl_permutation;

gsl_permutation *gsl_permutation_alloc (const size_t n);
gsl_permutation *gsl_permutation_calloc (const size_t n);
void gsl_permutation_init (gsl_permutation * p);
void gsl_permutation_free (gsl_permutation * p);
int gsl_permutation_memcpy (gsl_permutation * dest, const gsl_permutation * src);

int gsl_permutation_fread (FILE * istream, gsl_permutation * p);
int gsl_permutation_fwrite (FILE * ostream, const gsl_permutation * p);
int gsl_permutation_fscanf (FILE * istream, gsl_permutation * p);
int gsl_permutation_fprintf (FILE * ostream, const gsl_permutation * p, const char *format);

size_t gsl_permutation_size (const gsl_permutation * p);
size_t * gsl_permutation_data (const gsl_permutation * p);

int gsl_permutation_swap (gsl_permutation * p, const size_t i, const size_t j);

int gsl_permutation_valid (const gsl_permutation * p);
void gsl_permutation_reverse (gsl_permutation * p);
int gsl_permutation_inverse (gsl_permutation * inv, const gsl_permutation * p);
int gsl_permutation_next (gsl_permutation * p);
int gsl_permutation_prev (gsl_permutation * p);
int gsl_permutation_mul (gsl_permutation * p, const gsl_permutation * pa, const gsl_permutation * pb);

int gsl_permutation_linear_to_canonical (gsl_permutation * q, const gsl_permutation * p);
int gsl_permutation_canonical_to_linear (gsl_permutation * p, const gsl_permutation * q);

size_t gsl_permutation_inversions (const gsl_permutation * p);
size_t gsl_permutation_linear_cycles (const gsl_permutation * p);
size_t gsl_permutation_canonical_cycles (const gsl_permutation * q);

size_t gsl_permutation_get (const gsl_permutation * p, const size_t i);
