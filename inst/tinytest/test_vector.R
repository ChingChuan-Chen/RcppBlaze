## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

cppFile <- "test-vectors.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  Rcpp::sourceCpp(file.path("cpp", cppFile))
} else {
  Rcpp::sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
}
library(Matrix)
library(MatrixExtra)

vector_wrap_res <- vector_wrap_test()
expect_dbl_vec <- c(1.5, -2.5, 4.5)
expect_int_sparse_vector <- Matrix(c(0, 0, 1, 0, 3, 0), 6, 1)
expect_dbl_sparse_vector <- Matrix(c(0, 0, 1.5, 0, 3.6, 0), 6, 1)
expect_equal(vector_wrap_res[["dv_int"]], c(1L, -2L, 4L), info = "dv_int")
expect_equal(vector_wrap_res[["dv_cmpl"]], c(1.5+0.2*1i, -0.75-0.3*1i, -3+0.1*1i), info = "dv_cmpl")
expect_equal(vector_wrap_res[["dv_double"]], expect_dbl_vec, info = "dv_double")
expect_equal(vector_wrap_res[["sv_double"]], expect_dbl_vec, info = "sv_double")
expect_equal(vector_wrap_res[["hv_double"]], expect_dbl_vec, info = "hv_double")
expect_equal(vector_wrap_res[["cv_ua_up_double"]], expect_dbl_vec, info = "cv_ua_up_double")
expect_equal(vector_wrap_res[["cv_ua_pa_double"]], expect_dbl_vec, info = "cv_ua_pa_double")
expect_equal(vector_wrap_res[["cv_al_up_double"]], expect_dbl_vec, info = "cv_al_up_double")
expect_equal(vector_wrap_res[["cv_al_pa_double"]], expect_dbl_vec, info = "cv_al_pa_double")
expect_equal(vector_wrap_res[["uv_double"]], rep(3, 6L), info = "uv_double")
expect_equal(vector_wrap_res[["uv_cplx"]], rep(3.3-2.7*1i, 6L), info = "uv_cplx")
expect_equal(vector_wrap_res[["cpv_int"]], expect_int_sparse_vector, info = "cpv_int")
expect_equal(vector_wrap_res[["cpv_double"]], expect_dbl_sparse_vector, info = "cpv_double")
expect_equal(vector_wrap_res[["cpv_double_rv"]], as.csc.matrix(t(expect_dbl_sparse_vector)), info = "cpv_double_rv")
expect_equal(vector_wrap_res[["cpv_float"]], expect_dbl_sparse_vector, info = "cpv_float", tolerance = 1e-6)
expect_equal(vector_wrap_res[["zv_double"]], Matrix(0, 6L, 1L), info = "zv_double")

vector_as_res <- vector_as_test(list(c(1L, 3L, 6L), c(1.5, 2.5, 4.5)))
expect_double_sum <- 8.5
expect_equal(vector_as_res[["dv_int_sum"]], 10L, info = "dv_int_sum")
expect_equal(vector_as_res[["dv_double_sum"]], expect_double_sum, info = "dv_double_sum")
expect_equal(vector_as_res[["sv_double_sum"]], expect_double_sum, info = "sv_double_sum")
expect_equal(vector_as_res[["sv_double_unaligned_sum"]], expect_double_sum, info = "sv_double_unaligned_sum")
expect_equal(vector_as_res[["hv_double_sum"]], expect_double_sum, info = "hv_double_sum")
expect_equal(vector_as_res[["hv_double_unaligned_sum"]], expect_double_sum, info = "hv_double_unaligned_sum")

expect_error(vector_sv_error(c(1.5, 2.5, 4.5, 5.5)))
expect_error(vector_hv_error(c(1.5, 2.5, 4.5, 5.5)))

# custom_vector_as_res <- custom_vector_as_test(list(c(1L, 3L, 6L), c(1.5, 2.5, 4.5)))
# expect_equal(custom_vector_as_res[["iCustomVectorUU"]], 10L, info = "iCustomVectorUU")
# expect_equal(custom_vector_as_res[["dCustomVectorUU"]], expect_double_sum, info = "dCustomVectorUU")
# expect_equal(custom_vector_as_res[["dCustomVectorUP"]], expect_double_sum, info = "dCustomVectorUP")
# expect_equal(custom_vector_as_res[["dCustomVectorAU"]], expect_double_sum, info = "dCustomVectorAU")
# expect_equal(custom_vector_as_res[["dCustomVectorAP"]], expect_double_sum, info = "dCustomVectorAP")

sparse_vector_as_res <- sparse_vector_as_test(
  list(
    expect_int_sparse_vector,
    as.coo.matrix(expect_int_sparse_vector),
    as.csc.matrix(expect_int_sparse_vector),
    expect_dbl_sparse_vector,
    as.coo.matrix(expect_dbl_sparse_vector),
    as.csc.matrix(expect_dbl_sparse_vector),
    t(expect_dbl_sparse_vector),
    as.coo.matrix(t(expect_dbl_sparse_vector)),
    as.csc.matrix(t(expect_dbl_sparse_vector))
  )
)
expect_int_sum <- 4
expect_double_sum <- 5.1
expect_equal(sparse_vector_as_res[["cpv_int_dgCMatrix"]], expect_int_sum, info = "cpv_int_dgCMatrix")
expect_equal(sparse_vector_as_res[["cpv_int_dgTMatrix"]], expect_int_sum, info = "cpv_int_dgTMatrix")
expect_equal(sparse_vector_as_res[["cpv_int_dgRMatrix"]], expect_int_sum, info = "cpv_int_dgRMatrix")
expect_equal(sparse_vector_as_res[["cpv_double_dgCMatrix"]], expect_double_sum, info = "cpv_double_dgCMatrix")
expect_equal(sparse_vector_as_res[["cpv_double_dgTMatrix"]], expect_double_sum, info = "cpv_double_dgTMatrix")
expect_equal(sparse_vector_as_res[["cpv_double_dgRMatrix"]], expect_double_sum, info = "cpv_double_dgRMatrix")
expect_equal(sparse_vector_as_res[["cpv_double_rv_dgCMatrix"]], expect_double_sum, info = "cpv_double_rv_dgCMatrix")
expect_equal(sparse_vector_as_res[["cpv_double_rv_dgTMatrix"]], expect_double_sum, info = "cpv_double_rv_dgTMatrix")
expect_equal(sparse_vector_as_res[["cpv_double_rv_dgRMatrix"]], expect_double_sum, info = "cpv_double_rv_dgRMatrix")

expect_error(vector_cpv_wrong_row_error(Matrix(c(0, 0, 1), 3L, 1L)))
expect_error(vector_cpv_wrong_column_error(Matrix(c(0, 0, 1), 1L, 3L)))
