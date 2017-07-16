// Copyright (C)  2017         Chingchuan Chen
// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is based on RcppArmadillo.cpp and RcppEigen.h from RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.
//
// RcppBlaze.cpp: Rcpp/Blaze glue
//
// Copyright (C)  2017 Chingchuan Chen
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppBlaze is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppBlaze.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppBlaze.h>
#include <blaze/system/Version.h>

// [[Rcpp::export]]
Rcpp::IntegerVector blaze_version(bool single) {

  if (single)
    return Rcpp::wrap( 10 * BLAZE_MAJOR_VERSION + BLAZE_MINOR_VERSION );

  return Rcpp::IntegerVector::create(Rcpp::_["major"] = BLAZE_MAJOR_VERSION,
                                     Rcpp::_["minor"] = BLAZE_MINOR_VERSION );
}

// [[Rcpp::export]]
bool Blaze_SSE() {
#if defined(__MMX__) || defined(__SSE__) || defined(__SSE2__) || defined(__SSSE3__) || defined(__SSE4_1__) || defined(__SSE4_2__)
  return true;
#else
  return false;
#endif
}

// [[Rcpp::export]]
bool Blaze_AVX() {
#if defined(__AVX__)
  return true;
#else
  return false;
#endif
}

// [[Rcpp::export]]
bool Blaze_AVX2() {
#if defined(__AVX2__)
  return true;
#else
  return false;
#endif
}

// [[Rcpp::export]]
bool Blaze_MIC() {
#if defined(__MIC__)
  return true;
#else
  return false;
#endif
}

// [[Rcpp::export]]
bool Blaze_FMA() {
#if defined(__FMA__)
  return true;
#else
  return false;
#endif
}
