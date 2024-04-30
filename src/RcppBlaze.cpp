// Copyright (C)  2017 - 2024  Ching-Chuan Chen
//
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.

#include <RcppBlaze.h>
#include <blaze/system/Version.h>

//' The version of Blaze used in RcppBlaze
//'
//' To return the version of Blaze used in RcppBlaze.
//'
//' @param single A logical value indicates which type to return. If TRUE, it returns an integer. If FALSE, it returns a named vector.
//' @return A number or a named vector to represent the version of \code{blaze} depending on the input, \code{single}.
//' @seealso Blaze header file \code{blaze/system/Version.h}.
//' @examples
//' blaze_version(FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector blaze_version(bool single) {
  if (single) {
    return Rcpp::wrap(10 * BLAZE_MAJOR_VERSION + BLAZE_MINOR_VERSION);
  }

  return Rcpp::IntegerVector::create(
    Rcpp::_["major"] = BLAZE_MAJOR_VERSION,
    Rcpp::_["minor"] = BLAZE_MINOR_VERSION
  );
}
