/** \file
\brief Compute ECDF of a set of points (from stdin), puts results on stdout
 */

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/**Count number of file lines*/
int file_nlines(char* fname) {
  FILE* file = fopen(fname, "r");
  int ch, prev = '\n' /* so empty files have no lines */, lines = 0;
  while ( (ch = fgetc(file)) != EOF ) {/* Read all chars in the file. */
    if ( ch == '\n' ){
      ++lines; /* Bump the counter for every newline. */
    }
    prev = ch; /* Keep a copy to later test whether... */
  }
  fclose(file);
  if ( prev != '\n' ) /* ...the last line did not end in a newline. */
      ++lines; /* If so, add one more to the total. */
  return lines;
}


/**Update vector of Fn values recursively*/
void Fn_update(gsl_vector* Fn, gsl_vector* xn, int n, gsl_matrix* X0,
	       gsl_vector* workspace) {
  gsl_vector_scale(Fn, (double) n);
  for(int i=0; i < X0->size1; i++) {
    gsl_vector_memcpy(workspace, xn);
    gsl_vector_view zv = gsl_matrix_row(X0, i);
    gsl_vector_sub(workspace, &(zv.vector));
    gsl_vector_set(Fn, i, gsl_vector_get(Fn, i) + (!gsl_vector_ispos(workspace)));
  }
  gsl_vector_scale(Fn, 1.0 / ((double) n + 1.0));
}

int main(int argc, char** argv) {
  if(argc != 3) {
    printf("Compute ECDF from a sample\n");
    printf("usage: %s dim X0_filename\n", argv[0]);
    printf("\tdim: state space dimension\n");
    printf("\tX0_filename: file containing points in which to compute the ECDF\n");
    printf("reads input from stdin, puts output on stdout\n");
    return 1;
  }
  gsl_set_error_handler_off();
  int dim;
  sscanf(argv[1], "%d", &dim);
  int X0_N = file_nlines(argv[2]) / dim;
  FILE* f = fopen(argv[2], "r");
  gsl_matrix* X0 = gsl_matrix_alloc(X0_N, dim);
  gsl_matrix_fscanf(f, X0);
  fclose(f);
  gsl_vector* Fn = gsl_vector_alloc(X0_N);
  gsl_vector_set_all(Fn, 0.0);
  gsl_vector* xn = gsl_vector_alloc(dim);
  gsl_vector* workspace = gsl_vector_alloc(dim);
  int n=0;
  while(!gsl_vector_fscanf(stdin, xn)) {
    Fn_update(Fn, xn, n, X0, workspace);
    n++;
  }
  gsl_vector_fprintf(stdout, Fn, "%f");
  gsl_vector_free(workspace);
  gsl_vector_free(xn);
  gsl_vector_free(Fn);
  gsl_matrix_free(X0);
  return 0;
}
