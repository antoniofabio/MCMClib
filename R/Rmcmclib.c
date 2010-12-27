#include <R.h>
#include <Rinternals.h>
#include <mcmclib/mh.h>

SEXP R_mcmclib_distribfun_alloc(SEXP f, SEXP rho) {
  SEXP ans;
  PROTECT(ans = allocVector(LISTSXP, 2));
  SET_VECTOR_ELT(ans, 0, f);
  SET_VECTOR_ELT(ans, 1, rho);
  UNPROTECT(1);
  return(ans);
}

double R_mcmclib_distribfun(void* data, const gsl_vector* x) {
  Rprintf("R_mcmclib_distrifun: %p\n", data);
  SEXP p = (SEXP) data;
  SEXP R_fcall, fn, rho, R_x;
  fn = VECTOR_ELT(p, 0);
  rho = VECTOR_ELT(p, 1);
  PROTECT(R_x = R_MakeExternalPtr(x, NULL, NULL));
  PROTECT(R_fcall = lang2(fn, R_NilValue));
  SETCADR(R_fcall, R_x);
  double ans = REAL(eval(R_fcall, rho))[0];
  UNPROTECT(2);
  return(ans);
}
