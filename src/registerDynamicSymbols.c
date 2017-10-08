#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _RcppBlaze_fastLmPure(SEXP, SEXP, SEXP);
extern SEXP _RcppBlaze_blaze_version(SEXP);
extern SEXP _RcppBlaze_Blaze_SSE();
extern SEXP _RcppBlaze_Blaze_AVX();
extern SEXP _RcppBlaze_Blaze_AVX2();
extern SEXP _RcppBlaze_Blaze_MIC();
extern SEXP _RcppBlaze_Blaze_FMA();

static const R_CallMethodDef CallEntries[] = {
  {"_RcppBlaze_fastLmPure", (DL_FUNC) &_RcppBlaze_fastLmPure, 3},
  {"_RcppBlaze_blaze_version", (DL_FUNC) &_RcppBlaze_blaze_version, 1},
  {"_RcppBlaze_Blaze_SSE", (DL_FUNC) &_RcppBlaze_Blaze_SSE, 0},
  {"_RcppBlaze_Blaze_AVX", (DL_FUNC) &_RcppBlaze_Blaze_AVX, 0},
  {"_RcppBlaze_Blaze_AVX2", (DL_FUNC) &_RcppBlaze_Blaze_AVX2, 0},
  {"_RcppBlaze_Blaze_MIC", (DL_FUNC) &_RcppBlaze_Blaze_MIC, 0},
  {"_RcppBlaze_Blaze_FMA", (DL_FUNC) &_RcppBlaze_Blaze_FMA, 0},
  {NULL, NULL, 0}
};

void R_init_RcppBlaze(DllInfo* info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
