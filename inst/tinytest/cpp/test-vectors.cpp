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
using namespace std::complex_literals;

// [[Rcpp::export]]
Rcpp::List vector_wrap_test() {
  blaze::DynamicVector<int> dv_int(3UL);
  dv_int = {1, -2, 4};

  blaze::DynamicVector<std::complex<double>> dv_cmpl(3UL);
  dv_cmpl = {1.5 + 0.2i, -0.75 - 0.3i, -3.0 + 0.1i};

  blaze::DynamicVector<double> dv_double(3UL);
  dv_double = {1.5, -2.5, 4.5};

  blaze::StaticVector<double, 3UL> sv_double;
  sv_double = {1.5, -2.5, 4.5};

  blaze::HybridVector<double, 3UL> hv_double(3UL);
  hv_double = {1.5, -2.5, 4.5};

  std::unique_ptr<double[], blaze::ArrayDelete> cv_ua_up_double_mem(new double[3UL]);
  blaze::CustomVector<double, blaze::unaligned, blaze::unpadded> cv_ua_up_double(cv_ua_up_double_mem.get(), 3UL);
  cv_ua_up_double = {1.5, -2.5, 4.5};

  std::unique_ptr<double[], blaze::ArrayDelete> cv_ua_pa_double_mem(new double[4UL]);
  blaze::CustomVector<double, blaze::unaligned, blaze::padded> cv_ua_pa_double(cv_ua_pa_double_mem.get(), 3UL, 4UL);
  cv_ua_pa_double = {1.5, -2.5, 4.5};

  std::unique_ptr<double[], blaze::Deallocate> cv_al_up_double_mem(blaze::allocate<double>(3UL));
  blaze::CustomVector<double, blaze::aligned, blaze::unpadded> cv_al_up_double(cv_al_up_double_mem.get(), 3UL);
  cv_al_up_double = {1.5, -2.5, 4.5};

  std::unique_ptr<double[], blaze::Deallocate> cv_al_pa_double_mem(blaze::allocate<double>(4UL));
  blaze::CustomVector<double, blaze::aligned, blaze::padded> cv_al_pa_double(cv_al_pa_double_mem.get(), 3UL, 4UL);
  cv_al_pa_double = {1.5, -2.5, 4.5};

  // TODO: ZeroVector

  blaze::CompressedVector<int> cv_int(6UL);
  cv_int[2] = 1;
  cv_int[4] = 3;

  blaze::CompressedVector<double> cv_double(6UL);
  cv_double[2] = 1.5;
  cv_double[4] = 3.6;

  blaze::CompressedVector<double, blaze::rowVector> cv_double_rv(6UL);
  cv_double_rv[2] = 1.5;
  cv_double_rv[4] = 3.6;

  blaze::CompressedVector<float> cv_float(6UL);
  cv_float[2] = 1.5;
  cv_float[4] = 3.6;

  blaze::ZeroVector<float> zv_double(6UL);

  return Rcpp::List::create(
    _["dv_int"] = dv_int,
    _["dv_cmpl"] = dv_cmpl,
    _["dv_double"] = dv_double,
    _["sv_double"] = sv_double,
    _["hv_double"] = hv_double,
    _["cv_ua_up_double"] = cv_ua_up_double,
    _["cv_ua_pa_double"] = cv_ua_pa_double,
    _["cv_al_up_double"] = cv_al_up_double,
    _["cv_al_pa_double"] = cv_al_pa_double,
    _["cv_int"] = cv_int,
    _["cv_double"] = cv_double,
    _["cv_double_rv"] = cv_double_rv,
    _["cv_float"] = cv_float,
    _["zv_double"] = zv_double
  );
}

// [[Rcpp::export]]
Rcpp::List vector_as_test(Rcpp::List input_list) {
  blaze::DynamicVector<double> dv_int = input_list[0];
  blaze::DynamicVector<double> dv_double = input_list[1];
  blaze::StaticVector<double, 3> sv_double = input_list[1];
  blaze::HybridVector<double, 3> hv_double = input_list[1];
  blaze::StaticVector<double, 3, blaze::rowVector, blaze::unaligned> sv_double_unaligned = input_list[1];
  blaze::HybridVector<double, 3, blaze::rowVector, blaze::unaligned> hv_double_unaligned = input_list[1];

  return Rcpp::List::create(
    _["dv_int_sum"] = blaze::sum(dv_int),
    _["dv_double_sum"] = blaze::sum(dv_double),
    _["sv_double_sum"] = blaze::sum(sv_double),
    _["sv_double_unaligned_sum"] = blaze::sum(sv_double_unaligned),
    _["hv_double_sum"] = blaze::sum(dv_double),
    _["hv_double_unaligned_sum"] = blaze::sum(hv_double_unaligned)
  );
}

// [[Rcpp::export]]
void vector_sv_error(Rcpp::NumericVector x) {
  blaze::StaticVector<double, 3> z = Rcpp::as<blaze::StaticVector<double, 3>>(x);
}

// [[Rcpp::export]]
void vector_hv_error(Rcpp::NumericVector x) {
  blaze::HybridVector<double, 3> z = Rcpp::as<blaze::HybridVector<double, 3>>(x);
}

// [[Rcpp::export]]
Rcpp::List custom_vector_as_test(Rcpp::List input_list) {
  typedef typename blaze::CustomVector<int, blaze::unaligned, blaze::unpadded> iCustomVectorUU;
  typedef typename blaze::CustomVector<double, blaze::unaligned, blaze::unpadded> dCustomVectorUU;
  typedef typename blaze::CustomVector<double, blaze::unaligned, blaze::padded> dCustomVectorUP;
  typedef typename blaze::CustomVector<double, blaze::aligned, blaze::unpadded> dCustomVectorAU;
  typedef typename blaze::CustomVector<double, blaze::aligned, blaze::padded> dCustomVectorAP;

  iCustomVectorUU cv_int = input_list[0];
  dCustomVectorUU cv_uu_double = input_list[1];
  dCustomVectorUP cv_up_double = input_list[1];
  dCustomVectorAU cv_au_double = input_list[1];
  dCustomVectorAP cv_ap_double = input_list[1];

  return Rcpp::List::create(
    Rcpp::_["iCustomVectorUU"] = blaze::sum(cv_int),
    Rcpp::_["dCustomVectorUU"] = blaze::sum(cv_uu_double),
    Rcpp::_["dCustomVectorUP"] = blaze::sum(cv_up_double),
    Rcpp::_["dCustomVectorAU"] = blaze::sum(cv_au_double),
    Rcpp::_["dCustomVectorAP"] = blaze::sum(cv_ap_double)
  );
}

/*
 CompressedVector, ZeroVector
*/
