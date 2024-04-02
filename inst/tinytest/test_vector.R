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
