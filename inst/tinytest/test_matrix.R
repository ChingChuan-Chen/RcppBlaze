## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.
suppressPackageStartupMessages({
  require(Rcpp)
  require(RcppBlaze)
  require(Matrix)
  require(MatrixExtra)
  require(tinytest)
})

cppFile <- "test-matrices.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  sourceCpp(file.path("cpp", cppFile))
} else {
  sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
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

matrix_wrap_res2 <- matrix_wrap_test2()
expect_int_sparse_matrix <- Matrix(c(0, 0, 2, 1, 0), 3L, 5L)
expect_dbl_sparse_matrix_cm <- Matrix(c(0, 0, 2.5, 1.5, 0), 3L, 5L)
expect_equal(matrix_wrap_res2[["um_double_cm"]], matrix(-3.2, 2L, 3L), info = "um_double_cm")
expect_equal(matrix_wrap_res2[["um_double_rm"]], matrix(-3.2, 2L, 3L), info = "um_double_rm")
expect_equal(matrix_wrap_res2[["um_cplx"]], matrix(-1.8 + 0.6 * 1i, 2L, 3L), info = "um_cplx")
expect_equal(matrix_wrap_res2[["cpm_int"]], expect_int_sparse_matrix, info = "cpm_int")
expect_equal(matrix_wrap_res2[["cpm_double_cm"]], expect_dbl_sparse_matrix_cm, info = "cpm_double_cm")
expect_equal(matrix_wrap_res2[["cpm_double_rm"]], as.csr.matrix(expect_dbl_sparse_matrix_cm), info = "cpm_double_rm")
expect_equal(matrix_wrap_res2[["cpm_float"]], expect_dbl_sparse_matrix_cm, info = "cpm_float", tolerance = 1e-6)
expect_equal(matrix_wrap_res2[["im_double_cm"]], Diagonal(3L, rep(1, 3L)), info = "im_double_cm")
expect_equal(matrix_wrap_res2[["im_double_rm"]], Diagonal(3L, rep(1, 3L)), info = "im_double_rm")
expect_equal(matrix_wrap_res2[["zm_double_cm"]], Matrix(0, 2L, 3L), info = "zm_double_cm")
expect_equal(matrix_wrap_res2[["zm_double_rm"]], Matrix(0, 2L, 3L), info = "zm_double_rm")

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

expect_error(matrix_sm_error(matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)))
expect_error(matrix_hm_error(matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)))

matrix_as_res <- custom_matrix_as_test(
  list(
    matrix(c(1L, 2L, 4L, 5L, 8L, 3L), nrow=2L, byrow=TRUE),
    matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.3), nrow=2L, byrow=TRUE)
  )
)
expect_equal(matrix_as_res[["iCustomMatrixUU"]], 23L, info = "iCustomMatrixUU")
expect_equal(matrix_as_res[["iCustomMatrixAP"]], 23L, info = "iCustomMatrixAP")
expect_equal(matrix_as_res[["iCustomMatrixAP_RM"]], 23L, info = "iCustomMatrixAP_RM")
expect_equal(matrix_as_res[["dCustomMatrixUU"]], expect_double_sum, info = "dCustomMatrixUU")
expect_equal(matrix_as_res[["dCustomMatrixUP"]], expect_double_sum, info = "dCustomMatrixUP")
expect_equal(matrix_as_res[["dCustomMatrixAU"]], expect_double_sum, info = "dCustomMatrixAU")
expect_equal(matrix_as_res[["dCustomMatrixAP"]], expect_double_sum, info = "dCustomMatrixAP")

# Column-Major SparseMatrix
expect_double_sm_sum <- 12
dgCMatrix_as_test_res <- sparse_matrix_as_test(list(expect_dbl_sparse_matrix_cm))
expect_equal(dgCMatrix_as_test_res[["cpm_cm"]], 12, info="dgCMatrix_cpm_cm")
expect_equal(dgCMatrix_as_test_res[["cpm_rm"]], 12, info="dgCMatrix_cpm_rm")

dims <- c(3L, 3L)
value <- c(3.25, 1.25, 1.5)
idx_i <-c(0L,0L,1L)
idx_j <- c(1L,2L,2L)
cm_p_vec <- c(0L,0L,1L,3L)
rm_p_vec <- c(0L,2L,3L,3L)
tri_sm_u <- new("dtCMatrix", x=value, i=idx_i, p=cm_p_vec, Dim=dims, uplo="U")
tri_sm_l <- new("dtCMatrix", x=value, i=idx_j, p=rm_p_vec, Dim=dims, uplo="L")
dtCMatrix_upper_as_test_res <- sparse_matrix_as_test(list(tri_sm_u))
dtCMatrix_lower_as_test_res <- sparse_matrix_as_test(list(tri_sm_l))
expect_equal(dtCMatrix_upper_as_test_res[["cpm_cm"]], 6, info="dtCMatrix_upper_cpm_cm")
expect_equal(dtCMatrix_upper_as_test_res[["cpm_rm"]], 6, info="dtCMatrix_upper_cpm_rm")
expect_equal(dtCMatrix_lower_as_test_res[["cpm_cm"]], 6, info="dtCMatrix_lower_cpm_cm")
expect_equal(dtCMatrix_lower_as_test_res[["cpm_rm"]], 6, info="dtCMatrix_lower_cpm_rm")

tri_sm_u_du <- new("dtCMatrix", x=value, i=idx_i, p=cm_p_vec, Dim=dims, uplo="U", diag="U")
tri_sm_l_du <- new("dtCMatrix", x=value, i=idx_j, p=rm_p_vec, Dim=dims, uplo="L", diag="U")
dtCMatrix_upper_unit_diag_as_test_res <- sparse_matrix_as_test(list(tri_sm_u_du))
dtCMatrix_lower_unit_diag_as_test_res <- sparse_matrix_as_test(list(tri_sm_l_du))
expect_equal(dtCMatrix_upper_unit_diag_as_test_res[["cpm_cm"]], 9, info="dtCMatrix_upper_unit_diag_cpm_cm")
expect_equal(dtCMatrix_upper_unit_diag_as_test_res[["cpm_rm"]], 9, info="dtCMatrix_upper_unit_diag_cpm_rm")
expect_equal(dtCMatrix_lower_unit_diag_as_test_res[["cpm_cm"]], 9, info="dtCMatrix_lower_unit_diag_cpm_cm")
expect_equal(dtCMatrix_lower_unit_diag_as_test_res[["cpm_rm"]], 9, info="dtCMatrix_lower_unit_diag_cpm_rm")

sym_sm_u <- new("dsCMatrix", x=value, i=idx_i, p=cm_p_vec, Dim=dims, uplo="U")
sym_sm_l <- new("dsCMatrix", x=value, i=idx_j, p=rm_p_vec, Dim=dims, uplo="L")
dsCMatrix_upper_as_test_res <- sparse_matrix_as_test(list(sym_sm_u))
dsCMatrix_lower_as_test_res <- sparse_matrix_as_test(list(sym_sm_l))
expect_equal(dsCMatrix_upper_as_test_res[["cpm_cm"]], 12, info="dsCMatrix_upper_cpm_cm")
expect_equal(dsCMatrix_upper_as_test_res[["cpm_rm"]], 12, info="dsCMatrix_upper_cpm_rm")
expect_equal(dsCMatrix_lower_as_test_res[["cpm_cm"]], 12, info="dsCMatrix_lower_cpm_cm")
expect_equal(dsCMatrix_lower_as_test_res[["cpm_rm"]], 12, info="dsCMatrix_lower_cpm_rm")

# Row-Major SparseMatrix
expect_double_sm_sum <- 12
dgRMatrix_as_test_res <- sparse_matrix_as_test(list(as.csr.matrix(expect_dbl_sparse_matrix_cm)))
expect_equal(dgRMatrix_as_test_res[["cpm_cm"]], 12, info="dgRMatrix_cpm_cm")
expect_equal(dgRMatrix_as_test_res[["cpm_rm"]], 12, info="dgRMatrix_cpm_rm")

tri_sm_u <- new("dtRMatrix", x=value, j=idx_j, p=rm_p_vec, Dim=dims, uplo="U")
tri_sm_l <- new("dtRMatrix", x=value, j=idx_i, p=cm_p_vec, Dim=dims, uplo="L")
dtRMatrix_upper_as_test_res <- sparse_matrix_as_test(list(tri_sm_u))
dtRMatrix_lower_as_test_res <- sparse_matrix_as_test(list(tri_sm_l))
expect_equal(dtRMatrix_upper_as_test_res[["cpm_cm"]], 6, info="dtRMatrix_upper_cpm_cm")
expect_equal(dtRMatrix_upper_as_test_res[["cpm_rm"]], 6, info="dtRMatrix_upper_cpm_rm")
expect_equal(dtRMatrix_lower_as_test_res[["cpm_cm"]], 6, info="dtRMatrix_lower_cpm_cm")
expect_equal(dtRMatrix_lower_as_test_res[["cpm_rm"]], 6, info="dtRMatrix_lower_cpm_rm")

tri_sm_u_du <- new("dtRMatrix", x=value, j=idx_j, p=rm_p_vec, Dim=dims, uplo="U", diag="U")
tri_sm_l_du <- new("dtRMatrix", x=value, j=idx_i, p=cm_p_vec, Dim=dims, uplo="L", diag="U")
dtRMatrix_upper_unit_diag_as_test_res <- sparse_matrix_as_test(list(tri_sm_u_du))
dtRMatrix_lower_unit_diag_as_test_res <- sparse_matrix_as_test(list(tri_sm_l_du))
expect_equal(dtRMatrix_upper_unit_diag_as_test_res[["cpm_cm"]], 9, info="dtRMatrix_upper_unit_diag_cpm_cm")
expect_equal(dtRMatrix_upper_unit_diag_as_test_res[["cpm_rm"]], 9, info="dtRMatrix_upper_unit_diag_cpm_rm")
expect_equal(dtRMatrix_lower_unit_diag_as_test_res[["cpm_cm"]], 9, info="dtRMatrix_lower_unit_diag_cpm_cm")
expect_equal(dtRMatrix_lower_unit_diag_as_test_res[["cpm_rm"]], 9, info="dtRMatrix_lower_unit_diag_cpm_rm")

sym_sm_u <- new("dsRMatrix", x=value, j=idx_j, p=rm_p_vec, Dim=dims, uplo="U")
sym_sm_l <- new("dsRMatrix", x=value, j=idx_i, p=cm_p_vec, Dim=dims, uplo="L")
dsRMatrix_upper_as_test_res <- sparse_matrix_as_test(list(sym_sm_u))
dsRMatrix_lower_as_test_res <- sparse_matrix_as_test(list(sym_sm_l))
expect_equal(dsRMatrix_upper_as_test_res[["cpm_cm"]], 12, info="dsRMatrix_upper_cpm_cm")
expect_equal(dsRMatrix_upper_as_test_res[["cpm_rm"]], 12, info="dsRMatrix_upper_cpm_rm")
expect_equal(dsRMatrix_lower_as_test_res[["cpm_cm"]], 12, info="dsRMatrix_lower_cpm_cm")
expect_equal(dsRMatrix_lower_as_test_res[["cpm_rm"]], 12, info="dsRMatrix_lower_cpm_rm")

# Index-based SparseMatrix
expect_double_sm_sum <- 12
dgTMatrix_as_test_res <- sparse_matrix_as_test(list(as.coo.matrix(expect_dbl_sparse_matrix_cm)))
expect_equal(dgTMatrix_as_test_res[["cpm_cm"]], 12, info="dgTMatrix_cpm_cm")
expect_equal(dgTMatrix_as_test_res[["cpm_rm"]], 12, info="dgTMatrix_cpm_rm")

tri_sm_u <- new("dtTMatrix", x=value, i=idx_i, j=idx_j, Dim=dims, uplo="U")
tri_sm_l <- new("dtTMatrix", x=value, j=idx_i, i=idx_j, Dim=dims, uplo="L")
dtTMatrix_upper_as_test_res <- sparse_matrix_as_test(list(tri_sm_u))
dtTMatrix_lower_as_test_res <- sparse_matrix_as_test(list(tri_sm_l))
expect_equal(dtTMatrix_upper_as_test_res[["cpm_cm"]], 6, info="dtTMatrix_upper_cpm_cm")
expect_equal(dtTMatrix_upper_as_test_res[["cpm_rm"]], 6, info="dtTMatrix_upper_cpm_rm")
expect_equal(dtTMatrix_lower_as_test_res[["cpm_cm"]], 6, info="dtTMatrix_lower_cpm_cm")
expect_equal(dtTMatrix_lower_as_test_res[["cpm_rm"]], 6, info="dtTMatrix_lower_cpm_rm")

tri_sm_u_du <- new("dtTMatrix", x=value, i=idx_i, j=idx_j, Dim=dims, uplo="U", diag="U")
tri_sm_l_du <- new("dtTMatrix", x=value, j=idx_i, i=idx_j, Dim=dims, uplo="L", diag="U")
dtTMatrix_upper_unit_diag_as_test_res <- sparse_matrix_as_test(list(tri_sm_u_du))
dtTMatrix_lower_unit_diag_as_test_res <- sparse_matrix_as_test(list(tri_sm_l_du))
expect_equal(dtTMatrix_upper_unit_diag_as_test_res[["cpm_cm"]], 9, info="dtTMatrix_upper_unit_diag_cpm_cm")
expect_equal(dtTMatrix_upper_unit_diag_as_test_res[["cpm_rm"]], 9, info="dtTMatrix_upper_unit_diag_cpm_rm")
expect_equal(dtTMatrix_lower_unit_diag_as_test_res[["cpm_cm"]], 9, info="dtTMatrix_lower_unit_diag_cpm_cm")
expect_equal(dtTMatrix_lower_unit_diag_as_test_res[["cpm_rm"]], 9, info="dtTMatrix_lower_unit_diag_cpm_rm")

sym_sm_u <- new("dsTMatrix", x=value, i=idx_i, j=idx_j, Dim=dims, uplo="U")
sym_sm_l <- new("dsTMatrix", x=value, j=idx_i, i=idx_j, Dim=dims, uplo="L")
dsTMatrix_symm_upper_as_test_res <- sparse_matrix_as_test(list(sym_sm_u))
dsTMatrix_symm_lower_as_test_res <- sparse_matrix_as_test(list(sym_sm_l))
expect_equal(dsTMatrix_symm_upper_as_test_res[["cpm_cm"]], 12, info="dsTMatrix_symm_upper_cpm_cm")
expect_equal(dsTMatrix_symm_upper_as_test_res[["cpm_rm"]], 12, info="dsTMatrix_symm_upper_cpm_rm")
expect_equal(dsTMatrix_symm_lower_as_test_res[["cpm_cm"]], 12, info="dsTMatrix_symm_lower_cpm_cm")
expect_equal(dsTMatrix_symm_lower_as_test_res[["cpm_rm"]], 12, info="dsTMatrix_symm_lower_cpm_rm")

# ddiMatrix
sparseDiagMat1 <- Diagonal(x = c(10, 1, 2, 0, 4))
ddiMatrix_as_test_res <- sparse_matrix_as_test(list(sparseDiagMat1))
expect_equal(ddiMatrix_as_test_res[["cpm_cm"]], 17, info="ddiMatrix_cpm_cm")
expect_equal(ddiMatrix_as_test_res[["cpm_rm"]], 17, info="ddiMatrix_cpm_rm")

sparseDiagMat2 <- new("ddiMatrix", Dim=c(5L, 5L), diag="U")
ddiMatrix_unit_diag_as_test_res <- sparse_matrix_as_test(list(sparseDiagMat2))
expect_equal(ddiMatrix_unit_diag_as_test_res[["cpm_cm"]], 5, info="ddiMatrix_unit_diag_cpm_cm")
expect_equal(ddiMatrix_unit_diag_as_test_res[["cpm_rm"]], 5, info="ddiMatrix_unit_diag_cpm_rm")
