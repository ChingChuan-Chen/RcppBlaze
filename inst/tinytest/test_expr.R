## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

cppFile <- "test-expr.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  Rcpp::sourceCpp(file.path("cpp", cppFile))
} else {
  Rcpp::sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
}
library(Matrix)
library(MatrixExtra)

vec_expr_wrap_res <- wrap_vec_expr_test()
dense_vec <- c(1.5, -2.5, 4.5)
sparse_vec <- new("dgCMatrix", x=c(1.3), i=c(0L), p=c(0L, 1L), Dim=c(3L, 1L))
dense_vec_added <- dense_vec + as.vector(sparse_vec)
expect_equal(vec_expr_wrap_res[["dv plus scalar"]], dense_vec + 2.0, info = "dv plus scalar")
expect_equal(vec_expr_wrap_res[["dv minus scalar"]], dense_vec - 2.0, info = "dv minus scalar")
expect_equal(vec_expr_wrap_res[["dv multiply scalar"]], dense_vec * 2.0, info = "dv multiply scalar")
expect_equal(vec_expr_wrap_res[["dv divide scalar"]], dense_vec / 2.0, info = "dv divide scalar")
expect_equal(vec_expr_wrap_res[["sv multiply scalar"]], sparse_vec * 2.0, info = "sv multiply scalar")
expect_equal(vec_expr_wrap_res[["sv divide scalar"]], sparse_vec / 2.0, info = "sv divide scalar")
expect_equal(vec_expr_wrap_res[["dv + dv"]], dense_vec + dense_vec, info = "dv + dv")
expect_equal(vec_expr_wrap_res[["dv + sv"]], dense_vec_added, info = "dv + sv")
expect_equal(vec_expr_wrap_res[["sv + dv"]], dense_vec_added, info = "sv + dv")

mat_expr_wrap_res <- wrap_mat_expr_test()
dense_mat <- matrix(c(1.5, -2.5, 4.5, -5.5, 8.5, -7.3), nrow=2L, byrow=TRUE)
sparse_mat <- as.csr.matrix(Matrix(c(0, 2.5, 1.5, 0, 0, 3.5), 2L, 3L))
dense_mat_added <- dense_mat + as.vector(sparse_mat)
expect_equal(mat_expr_wrap_res[["dm plus scalar"]], dense_mat + 2.0, info = "dm plus scalar")
expect_equal(mat_expr_wrap_res[["dm minus scalar"]], dense_mat - 2.0, info = "dm minus scalar")
expect_equal(mat_expr_wrap_res[["dm multiply scalar"]], dense_mat * 2.0, info = "dm multiply scalar")
expect_equal(mat_expr_wrap_res[["dm divide scalar"]], dense_mat / 2.0, info = "dm divide scalar")
expect_equal(mat_expr_wrap_res[["sm multiply scalar"]], sparse_mat * 2.0, info = "sm multiply scalar")
expect_equal(mat_expr_wrap_res[["sm divide scalar"]], sparse_mat / 2.0, info = "sm divide scalar")
expect_equal(mat_expr_wrap_res[["dm + dm"]], dense_mat + dense_mat, info = "dm + dm")
expect_equal(mat_expr_wrap_res[["dm + sm"]], dense_mat_added, info = "dm + sm")
expect_equal(mat_expr_wrap_res[["sm + dm"]], dense_mat_added, info = "sm + dm")
