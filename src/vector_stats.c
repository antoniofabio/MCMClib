#include "vector_stats.h"

void mcmclib_matrix_colmeans(gsl_matrix* m, gsl_vector* out) {
	gsl_vector_view cv;
	gsl_vector* col;
	for(int i=0; i<m->size2; i++) {
		cv = gsl_matrix_column(m, i);
		col = &(cv.vector);
		gsl_vector_set(out, i, gsl_stats_mean(col->data, col->stride, col->size));
	}
}

void mcmclib_matrix_rowmeans(gsl_matrix* m, gsl_vector* out) {
	gsl_vector_view rv;
	gsl_vector* row;
	for(int i=0; i<m->size1; i++) {
		rv = gsl_matrix_row(m, i);
		row = &(rv.vector);
		gsl_vector_set(out, i, gsl_stats_mean(row->data, row->stride, row->size));
	}
}
