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

//' Set/Get the random number generator for blaze with given seed
//'
//' @param seed A positive integer to specify the seed value for the random number generator.
//' @return No return value.
//' @rdname blaze_seed
//' @export
// [[Rcpp::export]]
void blaze_set_seed(uint32_t seed) {
  blaze::setSeed(seed);
}

//' @rdname blaze_seed
//' @export
// [[Rcpp::export]]
uint32_t blaze_get_seed() {
  return blaze::getSeed();
}

//' Set/Get the Number of Threads used in blaze
//'
//' @param n The number of threads to set in blaze.
//' @return \code{blaze_get_threads} returns an integer and \code{blaze_set_threads} returns nothing.
//' @seealso blaze wiki: \url{https://bitbucket.org/blaze-lib/blaze/wiki/Shared\%20Memory\%20Parallelization}.
//' @rdname blaze_threads
//' @export
// [[Rcpp::export]]
void blaze_set_num_threads(size_t n) {
#if BLAZE_HPX_PARALLEL_MODE || BLAZE_OPENMP_PARALLEL_MODE || BLAZE_CPP_THREADS_PARALLEL_MODE || BLAZE_BOOST_THREADS_PARALLEL_MODE
  blaze::setNumThreads(n);
#else
  (void) n;
#endif
}

//' @rdname blaze_threads
//' @export
// [[Rcpp::export]]
size_t blaze_get_num_threads() {
#if BLAZE_HPX_PARALLEL_MODE || BLAZE_OPENMP_PARALLEL_MODE || BLAZE_CPP_THREADS_PARALLEL_MODE || BLAZE_BOOST_THREADS_PARALLEL_MODE
  return blaze::getNumThreads();
#else
  return 1;
#endif
}
