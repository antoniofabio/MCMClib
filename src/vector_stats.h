#ifndef __VECTOR_STATS_H__
#define __VECTOR_STATS_H__
#include <gsl/gsl_statistics_double.h>
#include "common.h"

/** column means
*/
void mcmclib_matrix_colmeans(gsl_matrix* m, gsl_vector* out);
/** row means
*/
void mcmclib_matrix_rowmeans(gsl_matrix* m, gsl_vector* out);

/** get variance/covariance matrix out of the 'vertical' matrix 'm'
*/
void mcmclib_matrix_covariance(gsl_matrix* m, gsl_matrix* out);

#endif
