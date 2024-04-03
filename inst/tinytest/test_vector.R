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

cppFile <- "test-vectors.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  Rcpp::sourceCpp(file.path("cpp", cppFile))
} else {
  Rcpp::sourceCpp(system.file("tinytest", "cpp", cppFile, package = "RcppBlaze"))
}

vector_wrap_res <- vector_wrap_test()
expect_dbl_vec <- c(1.5, 2.5, 4.5)
expect_equal(vector_wrap_res[["dv_int"]], c(1L, 2L, 4L), info = "dv_int")
expect_equal(vector_wrap_res[["dv_cmpl"]], c(1.5+0.2*1i, 0.75+0.3*1i, 3+0.1*1i), info = "dv_cmpl")
expect_equal(vector_wrap_res[["dv_double"]], expect_dbl_vec, info = "dv_double")
expect_equal(vector_wrap_res[["sv_double"]], expect_dbl_vec, info = "sv_double")
expect_equal(vector_wrap_res[["hv_double"]], expect_dbl_vec, info = "hv_double")
expect_equal(vector_wrap_res[["cv_ua_up_double"]], expect_dbl_vec, info = "cv_ua_up_double")
expect_equal(vector_wrap_res[["cv_ua_pa_double"]], expect_dbl_vec, info = "cv_ua_pa_double")
expect_equal(vector_wrap_res[["cv_al_up_double"]], expect_dbl_vec, info = "cv_al_up_double")
expect_equal(vector_wrap_res[["cv_al_pa_double"]], expect_dbl_vec, info = "cv_al_pa_double")

vector_as_res <- vector_as_test(list(c(1L, 3L, 6L), c(1.5, 2.5, 4.5)))
expect_double_sum <- 8.5
expect_equal(vector_as_res[["dv_int_sum"]], 10L, info = "dv_int_sum")
expect_equal(vector_as_res[["dv_double_sum"]], expect_double_sum, info = "dv_double_sum")
expect_equal(vector_as_res[["sv_double_sum"]], expect_double_sum, info = "sv_double_sum")
expect_equal(vector_as_res[["sv_double_aligned_sum"]], expect_double_sum, info = "sv_double_aligned_sum")
expect_equal(vector_as_res[["hv_double_sum"]], expect_double_sum, info = "hv_double_sum")
expect_equal(vector_as_res[["hv_double_aligned_sum"]], expect_double_sum, info = "hv_double_aligned_sum")

expect_error(vector_sv_error(c(1.5, 2.5, 4.5, 5.5)))
expect_error(vector_hv_error(c(1.5, 2.5, 4.5, 5.5)))

vector_as_res <- custom_vector_as_test(list(c(1L, 3L, 6L), c(1.5, 2.5, 4.5)))
expect_equal(vector_as_res[["iCustomVectorUU"]], 10L, info = "iCustomVectorUU")
expect_equal(vector_as_res[["dCustomVectorUU"]], expect_double_sum, info = "dCustomVectorUU")
expect_equal(vector_as_res[["dCustomVectorUP"]], expect_double_sum, info = "dCustomVectorUP")
# expect_equal(vector_as_res[["dCustomVectorAU"]], expect_double_sum, info = "dCustomVectorAU")
# expect_equal(vector_as_res[["dCustomVectorAP"]], expect_double_sum, info = "dCustomVectorAP")

