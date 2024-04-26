// Copyright (C) 2017 - 2024 Ching-Chuan Chen
//
// This file is part of RcppBlaze.
//
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.

// [[Rcpp::depends(RcppBlaze)]]
#include <RcppBlaze.h>
using Rcpp::_;

// [[Rcpp::export]]
Rcpp::List wrap_vec_expr_test() {
  blaze::StaticVector<double, 3UL> dv;
  dv = {1.5, -2.5, 4.5};

  blaze::CompressedVector<double> sv(3UL);
  sv[0] = 1.3;

  return Rcpp::List::create(
    _["dv plus scalar"] = dv + 2.0,
    _["dv minus scalar"] = dv - 2.0,
    _["dv multiply scalar"] = dv * 2.0,
    _["dv divide scalar"] = dv / 2.0,
    // _["sv plus scalar"] = sv + 2.0,  // no such operator
    // _["sv minus scalar"] = sv - 2.0, // no such operator
    _["sv multiply scalar"] = sv * 2.0,
    _["sv divide scalar"] = sv / 2.0,
    _["dv + dv"] = dv + dv,
    _["dv + sv"] = dv + sv,
    _["sv + dv"] = sv + dv
    // _["sv + sv"] = sv + sv // unable to convert
  );
}

// [[Rcpp::export]]
Rcpp::List wrap_mat_expr_test() {
  blaze::StaticMatrix<double, 2UL, 3UL> dm;
  dm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::CompressedMatrix<double> sm(2UL, 3UL);
  sm(0, 1) = 1.5;
  sm(1, 0) = 2.5;
  sm(1, 2) = 3.5;

  return Rcpp::List::create(
    _["dm plus scalar"] = dm + 2.0,
    _["dm minus scalar"] = dm - 2.0,
    _["dm multiply scalar"] = dm * 2.0,
    _["dm divide scalar"] = dm / 2.0,
    // _["sm plus scalar"] = sv + 2.0,  // no such operator
    // _["sm minus scalar"] = sv - 2.0, // no such operator
    _["sm multiply scalar"] = sm * 2.0,
    _["sm divide scalar"] = sm / 2.0,
    _["dm + dm"] = dm + dm,
    _["dm + sm"] = dm + sm,
    _["sm + dm"] = sm + dm
    // _["sv + sv"] = sv + sv // unable to convert
  );
}
