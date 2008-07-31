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

void mcmclib_matrix_covariance(gsl_matrix* m, gsl_matrix* out) {
	int d = m->size2;
	int n = m->size1;
	gsl_matrix* mean = gsl_matrix_alloc(1, d);
	gsl_matrix_view rv;
	gsl_matrix* row;

	gsl_vector_view mv = gsl_matrix_row(mean, 0);
	mcmclib_matrix_colmeans(m, &(mv.vector));

	gsl_matrix_set_zero(out);
	for(int i=0; i<n; i++) {
		rv = gsl_matrix_submatrix (m, i, 0, 1, d);
		row = &(rv.matrix);
		gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, row, row, 1.0, out);
	}

	gsl_blas_dgemm (CblasTrans, CblasNoTrans, (double) -n, mean, mean, 1.0, out);

	gsl_matrix_scale(out, 1.0 / (double) n);

	gsl_matrix_free(mean);
}
