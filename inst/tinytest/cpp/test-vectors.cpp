// Copyright (C) 2017 - 2024 Ching-Chuan Chen
//
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppBlaze is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppBlaze If not, see <http://www.gnu.org/licenses/>.

// [[Rcpp::depends(RcppBlaze)]]
#include <RcppBlaze.h>
using Rcpp::_;
using namespace std::complex_literals;

// [[Rcpp::export]]
Rcpp::List vector_wrap_test() {

  blaze::DynamicVector<int> dv_int(3);
  dv_int[0] = 1;
  dv_int[1] = 2;
  dv_int[2] = 4;

  blaze::DynamicVector<std::complex<double>> dv_cmpl(3);
  dv_cmpl[0] = 1.5 + 0.2i;
  dv_cmpl[1] = 0.75 + 0.3i;
  dv_cmpl[2] = 3.0 + 0.1i;

  blaze::DynamicVector<double> dv_double(3);
  dv_double[0] = 1.5;
  dv_double[1] = 2.5;
  dv_double[2] = 4.5;

  blaze::StaticVector<double, 3> sv_double;
  sv_double[0] = 1.5;
  sv_double[1] = 2.5;
  sv_double[2] = 4.5;

  blaze::HybridVector<double, 3> hv_double(3);
  hv_double[0] = 1.5;
  hv_double[1] = 2.5;
  hv_double[2] = 4.5;

  std::vector<double> vec(3UL);
  blaze::CustomVector<double, blaze::unaligned, blaze::unpadded> cv_ua_up_double(&vec[0], 3UL);
  cv_ua_up_double[0] = 1.5;
  cv_ua_up_double[1] = 2.5;
  cv_ua_up_double[2] = 4.5;

  std::unique_ptr<double[]> cv_ua_pa_double_mem(new double[4]);
  blaze::CustomVector<double, blaze::unaligned, blaze::padded> cv_ua_pa_double(cv_ua_pa_double_mem.get(), 3UL, 4UL);
  cv_ua_pa_double[0] = 1.5;
  cv_ua_pa_double[1] = 2.5;
  cv_ua_pa_double[2] = 4.5;


  std::unique_ptr<double[], blaze::Deallocate> cv_al_up_double_mem(blaze::allocate<double>(3UL));
  blaze::CustomVector<double, blaze::aligned, blaze::unpadded> cv_al_up_double(cv_al_up_double_mem.get(), 3UL);
  cv_al_up_double[0] = 1.5;
  cv_al_up_double[1] = 2.5;
  cv_al_up_double[2] = 4.5;

  std::unique_ptr<double[], blaze::Deallocate> cv_al_pa_double_mem(blaze::allocate<double>(4UL));
  blaze::CustomVector<double, blaze::aligned, blaze::padded> cv_al_pa_double(cv_al_pa_double_mem.get(), 3UL, 4UL);
  cv_al_pa_double[0] = 1.5;
  cv_al_pa_double[1] = 2.5;
  cv_al_pa_double[2] = 4.5;

  // TODO: CompressedVector, ZeroVector

  return Rcpp::List::create(
    _["dv_int"] = dv_int,
    _["dv_cmpl"] = dv_cmpl,
    _["dv_double"] = dv_double,
    _["sv_double"] = sv_double,
    _["hv_double"] = hv_double,
    _["cv_ua_up_double"] = cv_ua_up_double,
    _["cv_ua_pa_double"] = cv_ua_pa_double,
    _["cv_al_up_double"] = cv_al_up_double,
    _["cv_al_pa_double"] = cv_al_pa_double
  );
}

/*

// [[Rcpp::export]]
Rcpp::List wrap_() {
  using namespace std::complex_literals;

  blaze::DynamicVector<double> dv_double(3);
  dv_double[0] = 1.5;
  dv_double[1] = 2.5;
  dv_double[2] = 4.5;

  blaze::DynamicVector<double, blaze::columnVector> dv_double1(3);
  dv_double1[0] = 1.5;
  dv_double1[1] = 2.5;
  dv_double1[2] = 4.5;

  blaze::DynamicVector<double, blaze::rowVector> dv_double2(3);
  dv_double2[0] = 1.5;
  dv_double2[1] = 2.5;
  dv_double2[2] = 4.5;

  // blaze::DynamicVector<std::complex<double>> dv_cplx({3.0 + 0.1i, 1.5 + 0.2i, 0.75 + 0.3i});
  //   blaze::DynamicVector<int> dp_int({3, 1, 2});
// blaze::DynamicVector<long> dp_long({3L, 1L, 2L});

  return Rcpp::List::create(
    Rcpp::_["blaze::DynamicVector<double>"] = dv_double,
    Rcpp::_["blaze::DynamicVector<double, blaze::columnVector>"] = dv_double1,
    Rcpp::_["blaze::DynamicVector<double, blaze::rowVector>"] = dv_double2
// Rcpp::_["blaze::DynamicVector<std::complex<double>>"] = dv_cplx,
//     Rcpp::_["blaze::DynamicVector<int>"] = dp_int,
// Rcpp::_["blaze::DynamicVector<long>"] = dp_long
  );
}

*/


// [[Rcpp::export]]
Rcpp::List asDynamicVector(Rcpp::NumericVector x1, Rcpp::IntegerVector x2, Rcpp::ComplexVector x3) {
  blaze::DynamicVector<double> dv_double = Rcpp::as<blaze::DynamicVector<double>>(x1);
  blaze::DynamicVector<double, blaze::rowVector> dv_double_r = Rcpp::as<blaze::DynamicVector<double, blaze::rowVector>>(x1);
  blaze::DynamicVector<int> dp_int = Rcpp::as<blaze::DynamicVector<int>>(x2);
  blaze::DynamicVector<std::complex<double>> dv_cplx = Rcpp::as<blaze::DynamicVector<std::complex<double>>>(x3);

  blaze::HybridVector<double, 10, blaze::rowVector> hv_double = Rcpp::as<blaze::HybridVector<double, 10, blaze::rowVector>>(x1);

  // StaticVector with default padding and alignment
  blaze::StaticVector<double, 10> sv_double1 = Rcpp::as<blaze::StaticVector<double, 10>>(x1);

  // StaticVector with default padding and alignment
  blaze::StaticVector<double, 10, blaze::rowVector> sv_double = Rcpp::as<blaze::StaticVector<double, 10, blaze::rowVector>>(x1);

  blaze::StaticVector<double, 10, blaze::rowVector, blaze::unaligned, blaze::unpadded> sv_double2 = Rcpp::as<blaze::StaticVector<double, 10, blaze::rowVector, blaze::unaligned, blaze::unpadded>>(x1);
  // Rcpp::Rcout << dv_double_r + hv_double + sv_double + sv_double2 << std::endl;

  // blaze::CustomVector<double, blaze::unaligned, blaze::unpadded, blaze::rowVector> cv_double = Rcpp::as<blaze::CustomVector<double, blaze::unaligned, blaze::unpadded, blaze::rowVector>>(x1);
  // Rcpp::Rcout << cv_double << std::endl;

  return Rcpp::List::create(
    Rcpp::_["blaze::DynamicVector<double>"] = hv_double
  );
}

/*
// [[Rcpp::export]]
blaze::StaticVector< std::complex<double>, 3UL > test_StaticVector_cpl_len3( blaze::StaticVector< std::complex<double>, 3UL > x) {
  return x;
}

// [[Rcpp::export]]
blaze::HybridVector<double, 3UL> test_HybridVector_dbl_len3( blaze::HybridVector<double, 3UL> x ) {
  return x;
}

// [[Rcpp::export]]
blaze::HybridVector< std::complex<double>, 3UL > test_HybridVector_cpl_len3( blaze::HybridVector< std::complex<double>, 3UL > x) {
  return x;
}

// [[Rcpp::export]]
blaze::DynamicVector<double> test_DynamicVector_dbl( blaze::DynamicVector<double> x) {
  return x;
}

// [[Rcpp::export]]
blaze::DynamicVector< std::complex<double> > test_DynamicVector_cpl( blaze::DynamicVector< std::complex<double> > x) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomVector< double, blaze::unaligned, blaze::unpadded >
  test_CustomVector1_dbl( blaze::CustomVector< double, blaze::unaligned, blaze::unpadded > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::CustomVector< double, blaze::aligned, blaze::unpadded >
  test_CustomVector2_dbl( blaze::CustomVector< double, blaze::aligned, blaze::unpadded > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::CustomVector< double, blaze::unaligned, blaze::padded >
  test_CustomVector3_dbl( blaze::CustomVector< double, blaze::unaligned, blaze::padded > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::CustomVector< double, blaze::aligned, blaze::padded >
  test_CustomVector4_dbl( blaze::CustomVector< double, blaze::aligned, blaze::padded > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::CustomVector< std::complex<double>, blaze::unaligned, blaze::unpadded >
  test_CustomVector1_cpl( blaze::CustomVector< std::complex<double>, blaze::unaligned, blaze::unpadded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomVector< std::complex<double>, blaze::aligned, blaze::unpadded >
  test_CustomVector2_cpl( blaze::CustomVector< std::complex<double>, blaze::aligned, blaze::unpadded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomVector< std::complex<double>, blaze::unaligned, blaze::padded >
  test_CustomVector3_cpl( blaze::CustomVector< std::complex<double>, blaze::unaligned, blaze::padded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomVector< std::complex<double>, blaze::aligned, blaze::padded >
  test_CustomVector4_cpl( blaze::CustomVector< std::complex<double>, blaze::aligned, blaze::padded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CompressedVector<double> test_CompressedVector_dbl_col( blaze::CompressedVector<double> x ) {
  return x;
}

// [[Rcpp::export]]
blaze::CompressedVector<double,blaze::rowVector>
  test_CompressedVector_dbl_row( blaze::CompressedVector<double,blaze::rowVector> x ) {
  return x;
}
*/
