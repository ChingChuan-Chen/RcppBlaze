## Copyright (C)       2017 Chingchuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppBlaze is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppBlaze  If not, see <http://www.gnu.org/licenses/>.

LdFlags <- function() {}

#' @importFrom Rcpp Rcpp.plugin.maker
inlineCxxPlugin <-  function() {
  getSettings <- Rcpp.plugin.maker(
    include.before = "#include <RcppBlaze.h>",
    libs = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CXXFLAGS)",
    package = c("BH", "RcppBlaze")
  )
  settings <- getSettings()
  settings$env$PKG_CXXFLAGS <- "$(SHLIB_OPENMP_CXXFLAGS)"
  settings$env$PKG_LIBS <- paste(settings$env$PKG_LIBS, LdFlags())
  return(settings)
}
