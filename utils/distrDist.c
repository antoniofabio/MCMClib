/** \file
\brief Compute distance from a target distribution

Target distrib. given as a sample of obs from stdin, puts output to stdout
Handles multiple parrallel chains as input
 */

#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>
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

/**Squares all elements of 'x'*/
void vector_square(gsl_vector* x) {
  for(int i=0; i < x->size; i++) {
    double xi = gsl_vector_get(x, i);
    gsl_vector_set(x, i, xi * xi);
  }
}

int main(int argc, char** argv) {
  if(argc != 4) {
    printf("Compute distance between two ECDF\n");
    printf("usage: %s M Fn_filename X0_filename\n", argv[0]);
    printf("\tM: number of X chains\n");
    printf("\tFn_filename: ...\n");
    printf("\tX0_filename: ...\n");
    printf("reads input from stdin\n");
    return 1;
  }
  gsl_set_error_handler_off();
  int M;
  sscanf(argv[1], "%d", &M);
  int X0_N = file_nlines(argv[2]);
  int dim = file_nlines(argv[3]) / X0_N;
  gsl_matrix* X0 = gsl_matrix_alloc(X0_N, dim);
  fflush(stdout);
  FILE *f = fopen(argv[3], "r");
  if(!f) {
    fprintf(stderr, "cannot open file %s", argv[3]);
    exit(1);
  }
  gsl_matrix_fscanf(f, X0);
  fclose(f);

  gsl_vector* Fn_Y = gsl_vector_alloc(X0_N);
  f = fopen(argv[2], "r");
  gsl_vector_fscanf(f, Fn_Y);
  fclose(f);

  gsl_vector* Fn_diff = gsl_vector_alloc(X0_N);
  gsl_vector* Fn_X[M];
  for(int m=0; m<M; m++) {
    Fn_X[m] = gsl_vector_alloc(X0_N);
    gsl_vector_set_all(Fn_X[m], 0.0);
  }
  gsl_vector* workspace = gsl_vector_alloc(dim);
  gsl_vector* xn = gsl_vector_alloc(dim);
  int n=0;

  /*main loop*/
  while(!feof(stdin)){
    for(int m=0; m<M && !gsl_vector_fscanf(stdin, xn); m++) {
      Fn_update(Fn_X[m], xn, n, X0, workspace);
      gsl_vector_memcpy(Fn_diff, Fn_X[m]);
      gsl_vector_sub(Fn_diff, Fn_Y);
      vector_square(Fn_diff);
      printf("%f\n", gsl_stats_mean(Fn_diff->data, Fn_diff->stride, Fn_diff->size));
    }
    n++;
  }
  /*end main loop*/

  gsl_vector_free(Fn_diff);
  for(int m=0; m<M; m++)
    gsl_vector_free(Fn_X[m]);
  gsl_vector_free(Fn_Y);
  gsl_vector_free(xn);
  gsl_vector_free(workspace);
  gsl_matrix_free(X0);
}
