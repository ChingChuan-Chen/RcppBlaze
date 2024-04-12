## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

cppFile <- "test-matrices.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  Rcpp::sourceCpp(file.path("cpp", cppFile))
} else {
  Rcpp::sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
}

matrix_wrap_res <- matrix_wrap_test()
expect_int_mat <- matrix(c(1L, -2L, 4L, 5L, -8L, 0L), nrow=2L, byrow=TRUE)
expect_cplx_mat <- matrix(c(1.5+0.2*1i, 2.5+0.4*1i, 4.5-0.5*1i, -0.5-0.1*1i, -12.1+2*1i, 0.6*1i), nrow=2L, byrow=TRUE)
expect_dbl_mat <- matrix(c(1.5, -2.5, 4.5, -5.5, 8.5, -7.3), nrow=2L, byrow=TRUE)
expect_equal(matrix_wrap_res[["dm_int_cm"]], expect_int_mat, info = "dm_int_cm")
expect_equal(matrix_wrap_res[["dm_int_rm"]], expect_int_mat, info = "dm_int_rm")
expect_equal(matrix_wrap_res[["dm_cmplx_cm"]], expect_cplx_mat, info = "dm_cmplx_cm")
expect_equal(matrix_wrap_res[["dm_cmplx_rm"]], expect_cplx_mat, info = "dm_cmplx_rm")
expect_equal(matrix_wrap_res[["dm_double_cm"]], expect_dbl_mat, info = "dm_double_cm")
expect_equal(matrix_wrap_res[["dm_double_rm"]], expect_dbl_mat, info = "dm_double_rm")
expect_equal(matrix_wrap_res[["sm_double_cm"]], expect_dbl_mat, info = "sm_double_cm")
expect_equal(matrix_wrap_res[["sm_double_rm"]], expect_dbl_mat, info = "sm_double_rm")
expect_equal(matrix_wrap_res[["hm_double_cm"]], expect_dbl_mat, info = "hm_double_cm")
expect_equal(matrix_wrap_res[["hm_double_rm"]], expect_dbl_mat, info = "hm_double_rm")
expect_equal(matrix_wrap_res[["cm_ua_up_double_cm"]], expect_dbl_mat, info = "cm_ua_up_double_cm")
expect_equal(matrix_wrap_res[["cm_ua_up_double_rm"]], expect_dbl_mat, info = "cm_ua_up_double_rm")
expect_equal(matrix_wrap_res[["cm_ua_pa_double_cm"]], expect_dbl_mat, info = "cm_ua_pa_double_cm")
expect_equal(matrix_wrap_res[["cm_ua_pa_double_rm"]], expect_dbl_mat, info = "cm_ua_pa_double_rm")
expect_equal(matrix_wrap_res[["cm_al_up_double_cm"]], expect_dbl_mat, info = "cm_al_up_double_cm")
expect_equal(matrix_wrap_res[["cm_al_up_double_rm"]], expect_dbl_mat, info = "cm_al_up_double_rm")
expect_equal(matrix_wrap_res[["cm_al_pa_double_cm"]], expect_dbl_mat, info = "cm_al_pa_double_cm")
expect_equal(matrix_wrap_res[["cm_al_pa_double_rm"]], expect_dbl_mat, info = "cm_al_pa_double_rm")
# add UniformMatrix wrap test

matrix_as_res <- matrix_as_test(
  list(
    matrix(c(1L, 2L, 4L, 5L, 8L, 3L), nrow=2L, byrow=TRUE),
    matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)
  )
)
expect_double_sum <- 29.8
expect_equal(matrix_as_res[["dm_int_sum"]], 23L, info = "dm_int_sum")
expect_equal(matrix_as_res[["dm_double_sum"]], expect_double_sum, info = "dm_double_sum")
expect_equal(matrix_as_res[["sm_double_sum"]], expect_double_sum, info = "sm_double_sum")
expect_equal(matrix_as_res[["sm_double_unaligned_sum"]], expect_double_sum, info = "sm_double_unaligned_sum")
expect_equal(matrix_as_res[["hm_double_sum"]], expect_double_sum, info = "hm_double_sum")
expect_equal(matrix_as_res[["hm_double_unaligned_sum"]], expect_double_sum, info = "hm_double_unaligned_sum")
# add UniformMatrix as test

expect_error(matrix_sm_error(matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)))
expect_error(matrix_hm_error(matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)))

matrix_as_res <- custom_matrix_as_test(
  list(
    matrix(c(1L, 2L, 4L, 5L, 8L, 3L), nrow=2L, byrow=TRUE),
    matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)
  )
)
expect_equal(matrix_as_res[["iCustomMatrixUU"]], 23L, info = "iCustomMatrixUU")
expect_equal(matrix_as_res[["dCustomMatrixUU"]], expect_double_sum, info = "dCustomMatrixUU")
expect_equal(matrix_as_res[["dCustomMatrixUP"]], expect_double_sum, info = "dCustomMatrixUP")
expect_equal(matrix_as_res[["dCustomMatrixAU"]], expect_double_sum, info = "dCustomMatrixAU")
expect_equal(matrix_as_res[["dCustomMatrixAP"]], expect_double_sum, info = "dCustomMatrixAP")

