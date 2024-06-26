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
Rcpp::List wrap_vector_view_test() {
  blaze::DynamicVector<double> dv_dbl(8UL);
  dv_dbl = {1.5, -2.3, 4.0, 4.4, -12.1, 8.6, 7.7, -5.9};

  blaze::CompressedVector<double> cpv_double(8UL);
  cpv_double[2] = 1.5;
  cpv_double[4] = 3.6;
  cpv_double[6] = 5.2;

  return Rcpp::List::create(
    _["subvector_dv_function"] = blaze::subvector(dv_dbl, 2UL, 4UL),
    _["subvector_dv_template"] = blaze::subvector<2UL, 4UL>(dv_dbl),
    _["elements_dv_function"] = blaze::elements(dv_dbl, {0UL, 2UL, 5UL}),
    _["elements_dv_template"] = blaze::elements<0UL, 2UL, 5UL>(dv_dbl),
    _["subvector_sv_function"] = blaze::subvector(cpv_double, 2UL, 4UL),
    _["subvector_sv_template"] = blaze::subvector<2UL, 4UL>(cpv_double),
    _["elements_sv_function"] = blaze::elements(cpv_double, {0UL, 2UL, 5UL}),
    _["elements_sv_template"] = blaze::elements<0UL, 2UL, 5UL>(cpv_double)
  );
}

#define GET_STATICMATRIX                                       \
  blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor> D{ \
  {3.0, 0.0, 0.0}, {8.0, 0.0, 0.0}, {-2.0, -1.0, 4.0}};

#define GET_COMPRESSEDMATRIX                                       \
  blaze::CompressedMatrix<double, blaze::rowMajor> E(3UL, 3UL); \
  E.reserve(5); E.append(0, 0, 3.0); E.finalize(0);                \
  E.append(1, 0, 8.0); E.finalize(1); E.append(2, 0, -2.0);        \
  E.append(2, 1, -1.0); E.append(2, 2, 4.0); E.finalize(2);

// [[Rcpp::export]]
Rcpp::List wrap_row_test() {
  GET_STATICMATRIX;

  auto row0 = blaze::row(D, 0UL);
  auto row1 = blaze::row(D, 1UL);
  auto row1_2 = blaze::row<1UL>(D);
  auto row2 = blaze::row(D, 2UL);

  GET_COMPRESSEDMATRIX;

  auto row3 = blaze::row(E, 0UL);
  auto row4 = blaze::row(E, 1UL);
  auto row4_2 = blaze::row<1UL>(E);
  auto row5 = blaze::row(E, 2UL);

  return Rcpp::List::create(
    _["row0_dm"] = row0,
    _["row1_dm"] = row1,
    _["row2_dm"] = row2,
    _["row1_dm_template"] = row1_2,
    _["row0_sm"] = row3,
    _["row1_sm"] = row4,
    _["row2_sm"] = row5,
    _["row1_sm_template"] = row4_2
  );
}

// [[Rcpp::export]]
Rcpp::List wrap_column_test() {
  GET_STATICMATRIX;

  auto column0 = blaze::column(D, 0UL);
  auto column1 = blaze::column(D, 1UL);
  auto column1_2 = blaze::column<1UL>(D);
  auto column2 = blaze::column(D, 2UL);

  GET_COMPRESSEDMATRIX;

  auto column3 = blaze::column(E, 0UL);
  auto column4 = blaze::column(E, 1UL);
  auto column4_2 = blaze::column<1UL>(E);
  auto column5 = blaze::column(E, 2UL);

  return Rcpp::List::create(
    _["column0_dm"] = column0,
    _["column1_dm"] = column1,
    _["column2_dm"] = column2,
    _["column1_dm_template"] = column1_2,
    _["column0_sm"] = column3,
    _["column1_sm"] = column4,
    _["column2_sm"] = column5,
    _["column1_sm_template"] = column4_2
  );
}


// [[Rcpp::export]]
Rcpp::List wrap_band_test() {
  GET_STATICMATRIX;

  auto band0 = blaze::band(D, 0L);
  auto band1 = blaze::band(D, -1L);
  auto band1_2 = blaze::band<-1L>(D);
  auto band2 = blaze::band(D, -2L);

  GET_COMPRESSEDMATRIX;

  auto band3 = blaze::band(E, 0L);
  auto band4 = blaze::band(E, -1L);
  auto band4_2 = blaze::band<-1L>(E);
  auto band5 = blaze::band(E, -2L);

  return Rcpp::List::create(
    _["band0_dm"] = band0,
    _["band1_dm"] = band1,
    _["band2_dm"] = band2,
    _["band1_dm_template"] = band1_2,
    _["band0_sm"] = band3,
    _["band1_sm"] = band4,
    _["band2_sm"] = band5,
    _["band1_sm_template"] = band4_2
  );
}

// [[Rcpp::export]]
Rcpp::List wrap_submatrix_test() {
  GET_STATICMATRIX;

  auto sm1 = blaze::submatrix<1UL, 0UL, 2UL, 2UL>(D);
  auto sm2 = blaze::submatrix(D, 1UL, 0UL, 2UL, 2UL);

  GET_COMPRESSEDMATRIX;

  auto sm3 = blaze::submatrix<1UL, 0UL, 2UL, 2UL>(E);
  auto sm4 = blaze::submatrix(E, 1UL, 0UL, 2UL, 2UL);

  return Rcpp::List::create(
    _["submat_dm_template"] = sm1,
    _["submat_dm_function"] = sm2,
    _["submat_sm_template"] = sm3,
    _["submat_sm_function"] = sm4
  );
}

// [[Rcpp::export]]
Rcpp::List wrap_rows_test() {
  GET_STATICMATRIX;

  auto rs1 = blaze::rows<0UL, 2UL>(D);
  auto rs2 = blaze::rows(D, {0UL, 2UL});

  GET_COMPRESSEDMATRIX;

  auto rs3 = blaze::rows<0UL, 2UL>(E);
  auto rs4 = blaze::rows(E, {0UL, 2UL});

  return Rcpp::List::create(
    _["rows_dm_template"] = rs1,
    _["rows_dm_function"] = rs2,
    _["rows_sm_template"] = rs3,
    _["rows_sm_function"] = rs4
  );
}

// [[Rcpp::export]]
Rcpp::List wrap_columns_test() {
  GET_STATICMATRIX;

  auto cs1 = blaze::columns<0UL, 2UL>(D);
  auto cs2 = blaze::columns(D, {0UL, 2UL});

  GET_COMPRESSEDMATRIX;

  auto cs3 = blaze::columns<0UL, 2UL>(E);
  auto cs4 = blaze::columns(E, {0UL, 2UL});

  return Rcpp::List::create(
    _["columns_dm_template"] = cs1,
    _["columns_dm_function"] = cs2,
    _["columns_sm_template"] = cs3,
    _["columns_sm_function"] = cs4
  );
}

