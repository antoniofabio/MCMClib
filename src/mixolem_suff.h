#ifndef __MCMCLIB_MIXOLEM_SUFF_H__
#define __MCMCLIB_MIXOLEM_SUFF_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/** \addtogroup misc
@{
\addtogroup mixolem
@{
*/


/** mixture complete data sufficient statistics
\note internal data structure
*/
typedef struct {
  gsl_vector* delta;
  gsl_vector** delta_x;
  gsl_matrix** delta_xx;
} mcmclib_mixolem_suff;

/**alloc a new mixolem_suff struct, pointing to externally owned data*/
mcmclib_mixolem_suff* mcmclib_mixolem_suff_makeref(gsl_vector* delta,
						   gsl_vector** delta_x,
						   gsl_matrix** delta_xx);
mcmclib_mixolem_suff* mcmclib_mixolem_suff_alloc(int K, int dim);
void mcmclib_mixolem_suff_free(mcmclib_mixolem_suff* p);
void mcmclib_mixolem_suff_memcpy(mcmclib_mixolem_suff* dest, mcmclib_mixolem_suff* src);
void mcmclib_mixolem_suff_add(mcmclib_mixolem_suff* dest, mcmclib_mixolem_suff* src);
void mcmclib_mixolem_suff_scale(mcmclib_mixolem_suff* p, double alpha);

/**@}*/
#endif
