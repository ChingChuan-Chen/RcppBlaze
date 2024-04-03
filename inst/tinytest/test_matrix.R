## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppBlaze is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppBlaze  If not, see <http://www.gnu.org/licenses/>.

cppFile <- "test-matrices.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  Rcpp::sourceCpp(file.path("cpp", cppFile))
} else {
  Rcpp::sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
}

matrix_wrap_res <- matrix_wrap_test()
expect_int_mat <- matrix(c(1L, 2L, 4L, 5L, 8L, 7L), nrow=2L)
expect_cplx_mat <- matrix(c(1.5+0.2*1i, 0.75+0.3*1i, 3+0.1*1i, 2.5+0.4*1i, 7.5+0.5*1i, 0.3+0.6*1i), nrow=2L)
expect_dbl_mat <- matrix(c(1.5, 2.5, 4.5, 5.5, 8.5, 7.5), nrow=2L)
expect_equal(matrix_wrap_res[["dm_int"]], expect_int_mat, info = "dm_int")
expect_equal(matrix_wrap_res[["dm_cmpl"]], expect_cplx_mat, info = "dm_cmpl")
expect_equal(matrix_wrap_res[["dm_double"]], expect_dbl_mat, info = "dm_double")
expect_equal(matrix_wrap_res[["sm_double"]], expect_dbl_mat, info = "sm_double")
expect_equal(matrix_wrap_res[["hm_double"]], expect_dbl_mat, info = "hm_double")
# expect_equal(matrix_wrap_res[["cm_ua_up_double"]], expect_dbl_mat, info = "cm_ua_up_double")
# expect_equal(matrix_wrap_res[["cm_ua_pa_double"]], expect_dbl_mat, info = "cm_ua_pa_double")
# expect_equal(matrix_wrap_res[["cm_al_up_double"]], expect_dbl_mat, info = "cm_al_up_double")
# expect_equal(matrix_wrap_res[["cm_al_pa_double"]], expect_dbl_mat, info = "cm_al_pa_double")

# vector_as_res <- vector_as_test(list(c(1L, 3L, 6L), c(1.5, 2.5, 4.5)))
# expect_double_sum <- 8.5
# expect_equal(vector_as_res[["dv_int_sum"]], 10L, info = "dv_int_sum")
# expect_equal(vector_as_res[["dv_double_sum"]], expect_double_sum, info = "dv_double_sum")
# expect_equal(vector_as_res[["sv_double_sum"]], expect_double_sum, info = "sv_double_sum")
# expect_equal(vector_as_res[["sv_double_aligned_sum"]], expect_double_sum, info = "sv_double_aligned_sum")
# expect_equal(vector_as_res[["hv_double_sum"]], expect_double_sum, info = "hv_double_sum")
# expect_equal(vector_as_res[["hv_double_aligned_sum"]], expect_double_sum, info = "hv_double_aligned_sum")
#
# vector_as_res <- custom_vector_as_test(list(c(1L, 3L, 6L), c(1.5, 2.5, 4.5)))
# expect_equal(vector_as_res[["iCustomVectorUU"]], 10L, info = "iCustomVectorUU")
# expect_equal(vector_as_res[["dCustomVectorUU"]], expect_double_sum, info = "dCustomVectorUU")
# expect_equal(vector_as_res[["dCustomVectorUP"]], expect_double_sum, info = "dCustomVectorUP")
# expect_equal(vector_as_res[["dCustomVectorAU"]], expect_double_sum, info = "dCustomVectorAU")
# expect_equal(vector_as_res[["dCustomVectorAP"]], expect_double_sum, info = "dCustomVectorAP")

