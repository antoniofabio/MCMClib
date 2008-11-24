/** \file
\brief Compute distance from a constant vector

Reads target vector from stdin, puts output to stdout
 */

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

int main(int argc, char** argv) {
  if(argc<=2)
    return 1;
  gsl_set_error_handler_off();
  FILE* in = fopen(argv[1], "r");
  int dim = argc - 2;
  gsl_vector* theta = gsl_vector_alloc(dim);
  gsl_vector* xi = gsl_vector_alloc(dim);
  while(!gsl_vector_fscanf(in, xi)) {
    gsl_vector_sub(xi, theta);
    double dist = 0.0;
    for(int i=0; i<dim; i++)
      dist += pow(gsl_vector_get(xi, i), 2.0);
    printf("%f\n", sqrt(dist) / (double) dim);
    fflush(stdout);
  }
  gsl_vector_free(xi);
  gsl_vector_free(theta);
  fclose(in);
}
