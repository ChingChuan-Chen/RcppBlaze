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
Rcpp::List matrix_wrap_test() {
  blaze::DynamicMatrix<int, blaze::columnMajor> dm_int_cm(2UL, 3UL);
  dm_int_cm = { { 1, -2, 4 }, { 5, -8, 0 } };

  blaze::DynamicMatrix<int, blaze::rowMajor> dm_int_rm(2UL, 3UL);
  dm_int_rm = { { 1, -2, 4 }, { 5, -8, 0 } };

  blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> dm_cmplx_cm(2, 3);
  dm_cmplx_cm = { {  1.5 + 0.2i,  2.5 + 0.4i, 4.5 - 0.5i }, { -0.5 - 0.1i, -12.1 + 2i, 0.6i } };

  blaze::DynamicMatrix<std::complex<double>, blaze::rowMajor> dm_cmplx_rm(2, 3);
  dm_cmplx_rm = { {  1.5 + 0.2i,  2.5 + 0.4i, 4.5 - 0.5i }, { -0.5 - 0.1i, -12.1 + 2i, 0.6i } };

  blaze::DynamicMatrix<double, blaze::columnMajor> dm_double_cm(2UL, 3UL);
  dm_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::DynamicMatrix<double, blaze::rowMajor> dm_double_rm(2UL, 3UL);
  dm_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::StaticMatrix<double, 2UL, 3UL, blaze::columnMajor> sm_double_cm;
  sm_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::StaticMatrix<double, 2UL, 3UL, blaze::rowMajor> sm_double_rm;
  sm_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::HybridMatrix<double, 2UL, 3UL, blaze::columnMajor> hm_double_cm(2UL, 3UL);
  hm_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::HybridMatrix<double, 2UL, 3UL, blaze::rowMajor> hm_double_rm(2UL, 3UL);
  hm_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::ArrayDelete> ua_up_memory_cm(new double[6UL]);
  blaze::CustomMatrix<double, blaze::unaligned, blaze::unpadded, blaze::columnMajor> cm_ua_up_double_cm(ua_up_memory_cm.get(), 2UL, 3UL);
  cm_ua_up_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::ArrayDelete> ua_up_memory_rm(new double[6UL]);
  blaze::CustomMatrix<double, blaze::unaligned, blaze::unpadded, blaze::rowMajor> cm_ua_up_double_rm(ua_up_memory_rm.get(), 2UL, 3UL);
  cm_ua_up_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::ArrayDelete> ua_pa_memory_cm(new double[6UL]);
  blaze::CustomMatrix<double, blaze::unaligned, blaze::padded, blaze::columnMajor> cm_ua_pa_double_cm(ua_pa_memory_cm.get(), 2UL, 3UL, 2UL);
  cm_ua_pa_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::ArrayDelete> ua_pa_memory_rm(new double[8UL]);
  blaze::CustomMatrix<double, blaze::unaligned, blaze::padded, blaze::rowMajor> cm_ua_pa_double_rm(ua_pa_memory_rm.get(), 2UL, 3UL, 4UL);
  cm_ua_pa_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::Deallocate> al_up_memory_cm(blaze::allocate<double>(6UL));
  blaze::CustomMatrix<double, blaze::aligned, blaze::unpadded, blaze::columnMajor> cm_al_up_double_cm(al_up_memory_cm.get(), 2UL, 3UL, 2UL);
  cm_al_up_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::Deallocate> al_up_memory_rm(blaze::allocate<double>(8UL));
  blaze::CustomMatrix<double, blaze::aligned, blaze::unpadded, blaze::rowMajor> cm_al_up_double_rm(al_up_memory_rm.get(), 2UL, 3UL, 4UL);
  cm_al_up_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::Deallocate> al_pa_memory_cm(blaze::allocate<double>(6UL));
  blaze::CustomMatrix<double, blaze::aligned, blaze::padded, blaze::columnMajor> cm_al_pa_double_cm(al_pa_memory_cm.get(), 2UL, 3UL, 2UL);
  cm_al_pa_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::Deallocate> al_pa_memory_rm(blaze::allocate<double>(8UL));
  blaze::CustomMatrix<double, blaze::aligned, blaze::padded, blaze::rowMajor> cm_al_pa_double_rm(al_pa_memory_rm.get(), 2UL, 3UL, 4UL);
  cm_al_pa_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  return Rcpp::List::create(
    _["dm_int_cm"] = dm_int_cm,
    _["dm_int_rm"] = dm_int_rm,
    _["dm_cmplx_cm"] = dm_cmplx_cm,
    _["dm_cmplx_rm"] = dm_cmplx_rm,
    _["dm_double_cm"] = dm_double_cm,
    _["dm_double_rm"] = dm_double_rm,
    _["sm_double_cm"] = sm_double_cm,
    _["sm_double_rm"] = sm_double_rm,
    _["hm_double_cm"] = hm_double_cm,
    _["hm_double_rm"] = hm_double_rm,
    _["cm_ua_up_double_cm"] = cm_ua_up_double_cm,
    _["cm_ua_up_double_rm"] = cm_ua_up_double_rm,
    _["cm_ua_pa_double_cm"] = cm_ua_pa_double_cm,
    _["cm_ua_pa_double_rm"] = cm_ua_pa_double_rm,
    _["cm_al_up_double_cm"] = cm_al_up_double_cm,
    _["cm_al_up_double_rm"] = cm_al_up_double_rm,
    _["cm_al_pa_double_cm"] = cm_al_pa_double_cm,
    _["cm_al_pa_double_rm"] = cm_al_pa_double_rm
  );
}

// [[Rcpp::export]]
Rcpp::List matrix_wrap_test2() {
  blaze::UniformMatrix<double, blaze::columnMajor> um_double_cm(2UL, 3UL, -3.2);
  blaze::UniformMatrix<double, blaze::rowMajor> um_double_rm(2UL, 3UL, -3.2);
  blaze::UniformMatrix<std::complex<double>, blaze::columnMajor> um_cplx(2UL, 3UL, -1.8 + 0.6i);

  blaze::CompressedMatrix<int, blaze::columnMajor> cpm_int(3UL, 5UL);
  cpm_int(0, 1) = 1;
  cpm_int(0, 4) = 2;
  cpm_int(1, 2) = 2;
  cpm_int(1, 4) = 1;
  cpm_int(2, 0) = 2;
  cpm_int(2, 2) = 1;

  blaze::CompressedMatrix<double, blaze::columnMajor> cpm_double_cm(3UL, 5UL);
  cpm_double_cm(0, 1) = 1.5;
  cpm_double_cm(0, 4) = 2.5;
  cpm_double_cm(1, 2) = 2.5;
  cpm_double_cm(1, 4) = 1.5;
  cpm_double_cm(2, 0) = 2.5;
  cpm_double_cm(2, 2) = 1.5;

  blaze::CompressedMatrix<double, blaze::rowMajor> cpm_double_rm(3UL, 5UL);
  cpm_double_rm(0, 1) = 1.5;
  cpm_double_rm(0, 4) = 2.5;
  cpm_double_rm(1, 2) = 2.5;
  cpm_double_rm(1, 4) = 1.5;
  cpm_double_rm(2, 0) = 2.5;
  cpm_double_rm(2, 2) = 1.5;

  blaze::CompressedMatrix<float, blaze::columnMajor> cpm_float(3UL, 5UL);
  cpm_float(0, 1) = 1.5;
  cpm_float(0, 4) = 2.5;
  cpm_float(1, 2) = 2.5;
  cpm_float(1, 4) = 1.5;
  cpm_float(2, 0) = 2.5;
  cpm_float(2, 2) = 1.5;

  blaze::IdentityMatrix<double, blaze::columnMajor> im_double_cm(3UL);
  blaze::IdentityMatrix<double, blaze::rowMajor> im_double_rm(3UL);

  blaze::ZeroMatrix<double, blaze::columnMajor> zm_double_cm(2UL, 3UL);
  blaze::ZeroMatrix<double, blaze::rowMajor> zm_double_rm(2UL, 3UL);

  return Rcpp::List::create(
    _["um_double_cm"] = um_double_cm,
    _["um_double_rm"] = um_double_rm,
    _["um_cplx"] = um_cplx,
    _["cpm_int"] = cpm_int,
    _["cpm_double_cm"] = cpm_double_cm,
    _["cpm_double_rm"] = cpm_double_rm,
    _["cpm_float"] = cpm_float,
    _["im_double_cm"] = im_double_cm,
    _["im_double_rm"] = im_double_rm,
    _["zm_double_cm"] = zm_double_cm,
    _["zm_double_rm"] = zm_double_rm
  );
}

// [[Rcpp::export]]
Rcpp::List matrix_as_test(Rcpp::List input_list) {
  blaze::DynamicMatrix<int, blaze::columnMajor> dm_int = input_list[0];
  blaze::DynamicMatrix<double, blaze::columnMajor> dm_double = input_list[1];
  blaze::StaticMatrix<double, 2, 3, blaze::columnMajor> sm_double = input_list[1];
  blaze::HybridMatrix<double, 2, 3, blaze::columnMajor> hm_double = input_list[1];
  blaze::StaticMatrix<double, 2, 3, blaze::columnMajor, blaze::unaligned> sm_double_unaligned = input_list[1];
  blaze::HybridMatrix<double, 2, 3, blaze::columnMajor, blaze::unaligned> hm_double_unaligned = input_list[1];

  return Rcpp::List::create(
    _["dm_int_sum"] = blaze::sum(dm_int),
    _["dm_double_sum"] = blaze::sum(dm_double),
    _["sm_double_sum"] = blaze::sum(sm_double),
    _["sm_double_unaligned_sum"] = blaze::sum(sm_double_unaligned),
    _["hm_double_sum"] = blaze::sum(dm_double),
    _["hm_double_unaligned_sum"] = blaze::sum(hm_double_unaligned)
  );
}

// [[Rcpp::export]]
Rcpp::List custom_matrix_as_test(Rcpp::List input_list) {
  using iCustomMatrixUU = blaze::CustomMatrix<int, blaze::unaligned, blaze::unpadded, blaze::columnMajor>;
  using iCustomMatrixAP = blaze::CustomMatrix<int, blaze::aligned, blaze::padded, blaze::columnMajor>;
  using iCustomMatrixAP_RM = blaze::CustomMatrix<int, blaze::aligned, blaze::padded, blaze::rowMajor>;

  // initialize
  Rcpp::Shield<SEXP> intMatDimsSexp(Rf_getAttrib(input_list[0], R_DimSymbol));
  int* intMatDims = INTEGER(intMatDimsSexp);
  size_t m = (size_t) intMatDims[0], n = (size_t) intMatDims[1];

  // column-major parameters
  size_t intSimdSize = blaze::SIMDTrait<int>::size;
  size_t intMatPaddedRows = blaze::nextMultiple<size_t>(m, intSimdSize);

  // unaligned & unpadded column-major CustomMatrix
  std::unique_ptr<int[], blaze::ArrayDelete> data_unpadded(new int[m*n]);
  iCustomMatrixUU cm_uu_int(data_unpadded.get(), m, n);
  RcppBlaze::copyToCustomMatrix(input_list[0], cm_uu_int);

  // aligned & padded column-major CustomMatrix
  std::unique_ptr<int[], blaze::Deallocate> data_padded(blaze::allocate<int>(intMatPaddedRows * n));
  iCustomMatrixAP cm_ap_int(data_padded.get(), m, n, intMatPaddedRows);
  RcppBlaze::copyToCustomMatrix(input_list[0], cm_ap_int);

  // row-major parameters
  size_t intMatPaddedCols = blaze::nextMultiple<size_t>(n, intSimdSize);

  // aligned & padded row-major CustomMatrix
  std::unique_ptr<int[], blaze::Deallocate> data_rm_padded(blaze::allocate<int>(m * intMatPaddedCols));
  iCustomMatrixAP_RM cm_ap_rm_int(data_rm_padded.get(), m, n, intMatPaddedCols);
  RcppBlaze::copyToCustomMatrix(input_list[0], cm_ap_rm_int);


  using dCustomMatrixUU = blaze::CustomMatrix<double, blaze::unaligned, blaze::unpadded, blaze::columnMajor>;
  using dCustomMatrixUP = blaze::CustomMatrix<double, blaze::unaligned, blaze::padded, blaze::columnMajor>;
  using dCustomMatrixAU = blaze::CustomMatrix<double, blaze::aligned, blaze::unpadded, blaze::columnMajor>;
  using dCustomMatrixAP = blaze::CustomMatrix<double, blaze::aligned, blaze::padded, blaze::columnMajor>;

  // initialize
  Rcpp::Shield<SEXP> dblMatDimsSexp(Rf_getAttrib(input_list[1], R_DimSymbol));
  int* dblMatDims = INTEGER(dblMatDimsSexp);
  size_t m2 = (size_t) dblMatDims[0], n2 = (size_t) dblMatDims[1];

  // column-major parameters
  size_t dblSimdSize = blaze::SIMDTrait<double>::size;
  size_t dblMatPaddedRows = blaze::nextMultiple<size_t>(m2, dblSimdSize);

  // unaligned & unpadded column-major CustomMatrix
  std::unique_ptr<double[], blaze::ArrayDelete> data_uu_dbl(new double[m2*n2]);
  dCustomMatrixUU cm_uu_dbl(data_uu_dbl.get(), m2, n2);
  RcppBlaze::copyToCustomMatrix(input_list[1], cm_uu_dbl);

  // unaligned & padded column-major CustomMatrix
  std::unique_ptr<double[], blaze::ArrayDelete> data_up_dbl(new double[dblMatPaddedRows * n2]);
  dCustomMatrixUP cm_up_dbl(data_up_dbl.get(), m2, n2, dblMatPaddedRows);
  RcppBlaze::copyToCustomMatrix(input_list[1], cm_up_dbl);

  // aligned & unpadded column-major CustomMatrix
  std::unique_ptr<double[], blaze::Deallocate> data_au_dbl(blaze::allocate<double>(m2*n2));
  dCustomMatrixAU cm_au_dbl(data_au_dbl.get(), m2, n2);
  RcppBlaze::copyToCustomMatrix(input_list[1], cm_au_dbl);

  // aligned & padded column-major CustomMatrix
  std::unique_ptr<double[], blaze::Deallocate> data_ap_dbl(blaze::allocate<double>(dblMatPaddedRows * n2));
  dCustomMatrixAP cm_ap_dbl(data_ap_dbl.get(), m2, n2, dblMatPaddedRows);
  RcppBlaze::copyToCustomMatrix(input_list[1], cm_ap_dbl);

  return Rcpp::List::create(
    Rcpp::_["iCustomMatrixUU"] = blaze::sum(cm_uu_int),
    Rcpp::_["iCustomMatrixAP"] = blaze::sum(cm_ap_int),
    Rcpp::_["iCustomMatrixAP_RM"] = blaze::sum(cm_ap_rm_int),
    Rcpp::_["dCustomMatrixUU"] = blaze::sum(cm_uu_dbl),
    Rcpp::_["dCustomMatrixUP"] = blaze::sum(cm_up_dbl),
    Rcpp::_["dCustomMatrixAU"] = blaze::sum(cm_au_dbl),
    Rcpp::_["dCustomMatrixAP"] = blaze::sum(cm_ap_dbl)
  );
}

// [[Rcpp::export]]
Rcpp::List sparse_matrix_as_test(Rcpp::List input_list) {
  blaze::CompressedMatrix<double, blaze::columnMajor> cpm_cm = input_list[0];
  blaze::CompressedMatrix<double, blaze::rowMajor> cpm_rm = input_list[0];

  return Rcpp::List::create(
    Rcpp::_["cpm_cm"] = blaze::sum(cpm_cm),
    Rcpp::_["cpm_rm"] = blaze::sum(cpm_rm)
  );
}
