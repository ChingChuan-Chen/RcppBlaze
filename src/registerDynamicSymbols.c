// Copyright (C) 2017 - 2024 Ching-Chuan Chen
//
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _RcppBlaze_blaze_version(SEXP);
extern SEXP _RcppBlaze_fastLmPure(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_RcppBlaze_blaze_version", (DL_FUNC) &_RcppBlaze_blaze_version, 1},
  {"_RcppBlaze_fastLmPure", (DL_FUNC) &_RcppBlaze_fastLmPure, 3},
  {NULL, NULL, 0}
};

void R_init_RcppBlaze(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
