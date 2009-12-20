%module mcmclib

%include "carrays.i"

%scheme %{(load-extension "libguilemcmclib.so" "scm_init_swig_mcmclib_module")%}

%{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

static void guile_mcmclib_err_handler(const char * reason,
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
  gsl_set_error_handler(guile_mcmclib_err_handler);
%}

%typemap(in) gsl_vector** {
  size_t $1_size = scm_c_vector_length($input);
  $1 = malloc($1_size * sizeof(gsl_vector*));
  for(int i=0; i<$1_size; i++) {
    $1[i] = (gsl_vector*) SWIG_MustGetPtr(scm_c_vector_ref($input, i),
					  SWIGTYPE_p_gsl_vector, $argnum, 0);
  }
}
%typemap(freearg) gsl_vector** {
  free($1);
}

%typemap(in) gsl_matrix** {
  size_t $1_size = scm_c_vector_length($input);
  $1 = malloc($1_size * sizeof(gsl_matrix*));
  for(int i=0; i<$1_size; i++) {
    $1[i] = (gsl_matrix*) SWIG_MustGetPtr(scm_c_vector_ref($input, i),
					  SWIGTYPE_p_gsl_matrix, $argnum, 0);
  }
}
%typemap(freearg) gsl_matrix** {
  free($1);
}

%include "mh_q.i"
%include "mh.i"
%include "amh.i"
%include "monitor.i"
%include "distrfuns.i"
%include "mixem.i"
