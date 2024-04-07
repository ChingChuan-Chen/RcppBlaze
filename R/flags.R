## Copyright (C) 2010 - 2024 Dirk Eddelbuettel, Romain Francois and Douglas Bates
## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is based on flags.R and inline.R from RcppArmadillo.
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

LdFlags <- function() {}

#' @importFrom Rcpp Rcpp.plugin.maker
inlineCxxPlugin <-  function() {
  getSettings <- Rcpp.plugin.maker(
    include.before = "#include <RcppBlaze.h>",
    libs = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)",
    package = c("RcppBlaze")
  )
  settings <- getSettings()
  settings$env$PKG_LIBS <- paste(settings$env$PKG_LIBS, LdFlags())
  return(settings)
}
