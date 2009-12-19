%module gsl

%include "carrays.i"

%scheme %{(load-extension "libguilegsl.so" "scm_init_swig_gsl_module")%}

%{
#include <string.h>
#include <libguile.h>
#include <gsl/gsl_errno.h>

static void guile_gsl_err_handler(const char * reason,
				  const char * file,
				  int line,
				  int gsl_errno) {
  static char msg[2048];
  sprintf(msg, "%s:%d", file, line);
  scm_error_scm(scm_misc_error_key,
		scm_from_locale_string(msg),
		scm_from_locale_string(reason), SCM_BOOL_F, SCM_BOOL_F);
}
%}

%init %{
  gsl_set_error_handler(guile_gsl_err_handler);
%}

%typemap(in) FILE* {
  $1 = fdopen(scm_to_int(scm_fileno($input)), "rw");
}
%typemap(freearg) FILE* {
  fclose($1);
}
%typemap(in) (FILE* istream) {
  $1 = fdopen(scm_to_int(scm_fileno($input)), "r");
}
%typemap(in) (FILE* ostream) {
  $1 = fdopen(scm_to_int(scm_fileno($input)), "w");
}

%include "vector.i"
%include "matrix.i"
%include "blas.i"
%include "rng.i"
%include "qrng.i"
%include "sf.i"
%include "mode.i"
%include "permutation.i"
%include "linalg.i"
