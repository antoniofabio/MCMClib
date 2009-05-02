%module mcmclib

%include "typemaps.i"
%include "carrays.i"
%include "gsl.i"

%array_functions(double, doubleArray);
%array_functions(char *, stringArray);
%array_functions(int, intArray);

%{
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lpdf_iwishart.h>
%}

mcmclib_iwishart_lpdf* mcmclib_iwishart_lpdf_alloc(gsl_matrix* Psi, int m);
void mcmclib_iwishart_lpdf_free(mcmclib_iwishart_lpdf* p);
double mcmclib_iwishart_lpdf_compute(void* p, gsl_vector* x);
