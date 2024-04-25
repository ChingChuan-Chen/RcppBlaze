## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

cppFile <- "test-views.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  Rcpp::sourceCpp(file.path("cpp", cppFile))
} else {
  Rcpp::sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
}
library(Matrix)

wrap_vector_views_res <- wrap_vector_view_test()
expect_equal(wrap_vector_views_res[["subvector_dv_function"]], c(4.0, 4.4, -12.1, 8.6), info = "subvector_dv_function")
expect_equal(wrap_vector_views_res[["subvector_dv_template"]], c(4.0, 4.4, -12.1, 8.6), info = "subvector_dv_template")
expect_equal(wrap_vector_views_res[["elements_dv_function"]], c(1.5, 4.0, 8.6), info = "elements_dv_function")
expect_equal(wrap_vector_views_res[["elements_dv_template"]], c(1.5, 4.0, 8.6), info = "elements_dv_template")
expect_equal(wrap_vector_views_res[["subvector_sv_function"]], Matrix(c(1.5, 0, 3.6, 0), 4, 1, sparse=TRUE), info = "subvector_sv_function")
expect_equal(wrap_vector_views_res[["subvector_sv_template"]], Matrix(c(1.5, 0, 3.6, 0), 4, 1, sparse=TRUE), info = "subvector_sv_template")
expect_equal(wrap_vector_views_res[["elements_sv_function"]], Matrix(c(0, 1.5, 0), 3, 1, sparse=TRUE), info = "elements_sv_function")
expect_equal(wrap_vector_views_res[["elements_sv_template"]], Matrix(c(0, 1.5, 0), 3, 1, sparse=TRUE), info = "elements_sv_template")

wrap_band_res <- wrap_band_test()
band0 <- c(3, 0, 4)
band1 <- c(8, -1)
band2 <- c(-2)
expect_equal(wrap_band_res[["band0_dm"]], band0, info = "band0_dm")
expect_equal(wrap_band_res[["band1_dm"]], band1, info = "band1_dm")
expect_equal(wrap_band_res[["band2_dm"]], band2, info = "band2_dm")
expect_equal(wrap_band_res[["band1_dm_template"]], band1, info = "band1_dm_template")
expect_equal(wrap_band_res[["band0_sm"]], Matrix(band0, 3, 1, sparse=TRUE), info = "band0_sm")
expect_equal(wrap_band_res[["band1_sm"]], Matrix(band1, 2, 1, sparse=TRUE), info = "band1_sm")
expect_equal(wrap_band_res[["band2_sm"]], new("dgCMatrix", x=band2, i=c(0L), p=c(0L, 1L), Dim=c(1L, 1L)), info = "band2_sm")
expect_equal(wrap_band_res[["band1_sm_template"]], Matrix(band1, 2, 1, sparse=TRUE), info = "band1_sm_template")

wrap_row_res <- wrap_row_test()
row0 <- c(3, 0, 0)
row1 <- c(8, 0, 0)
row2 <- c(-2, -1, 4)
expect_equal(wrap_row_res[["row0_dm"]], row0, info = "row0_dm")
expect_equal(wrap_row_res[["row1_dm"]], row1, info = "row1_dm")
expect_equal(wrap_row_res[["row2_dm"]], row2, info = "row2_dm")
expect_equal(wrap_row_res[["row1_dm_template"]], row1, info = "row1_dm_template")
expect_equal(wrap_row_res[["row0_sm"]], Matrix(row0, 1, 3, sparse=TRUE), info = "row0_sm")
expect_equal(wrap_row_res[["row1_sm"]], Matrix(row1, 1, 3, sparse=TRUE), info = "row1_sm")
expect_equal(wrap_row_res[["row2_sm"]], Matrix(row2, 1, 3, sparse=TRUE), info = "row2_sm")
expect_equal(wrap_row_res[["row1_sm_template"]], Matrix(row1, 1, 3, sparse=TRUE), info = "row1_sm_template")

wrap_column_res <- wrap_column_test()
column0 <- c(3, 8, -2)
column1 <- c(0, 0, -1)
column2 <- c(0, 0, 4)
expect_equal(wrap_column_res[["column0_dm"]], column0, info = "column0_dm")
expect_equal(wrap_column_res[["column1_dm"]], column1, info = "column1_dm")
expect_equal(wrap_column_res[["column2_dm"]], column2, info = "column2_dm")
expect_equal(wrap_column_res[["column1_dm_template"]], column1, info = "column1_dm_template")
expect_equal(wrap_column_res[["column0_sm"]], Matrix(column0, 3, 1, sparse=TRUE), info = "column0_sm")
expect_equal(wrap_column_res[["column1_sm"]], Matrix(column1, 3, 1, sparse=TRUE), info = "column1_sm")
expect_equal(wrap_column_res[["column2_sm"]], Matrix(column2, 3, 1, sparse=TRUE), info = "column2_sm")
expect_equal(wrap_column_res[["column1_sm_template"]], Matrix(column1, 3, 1, sparse=TRUE), info = "column1_sm_template")
