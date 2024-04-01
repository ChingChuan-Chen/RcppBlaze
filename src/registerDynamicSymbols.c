#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _RcppBlaze_blaze_version(SEXP);
extern SEXP _RcppBlaze_testAs1(SEXP);
extern SEXP _RcppBlaze_testAs2(SEXP);
extern SEXP _RcppBlaze_testAs3(SEXP);
extern SEXP _RcppBlaze_testAs4(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_RcppBlaze_blaze_version", (DL_FUNC) &_RcppBlaze_blaze_version, 1},
  {"_RcppBlaze_testAs1", (DL_FUNC) &_RcppBlaze_testAs1, 1},
  {"_RcppBlaze_testAs2", (DL_FUNC) &_RcppBlaze_testAs2, 1},
  {"_RcppBlaze_testAs3", (DL_FUNC) &_RcppBlaze_testAs3, 1},
  {"_RcppBlaze_testAs4", (DL_FUNC) &_RcppBlaze_testAs4, 1},
  {NULL, NULL, 0}
};

void R_init_RcppBlaze(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
