## Copyright (C)       2017 Chingchuan Chen
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
  sourceCpp(file.path( "cpp", cppFile))
} else {
  sourceCpp(system.file("unitTests", "cpp", cppFile, package = "RcppBlaze"))
}

testVectors <- function(){
  x <- 1:3
  z <- 1:3 + 2i
  sv_col <- Matrix::sparseMatrix(i = c(3, 4), j = c(1, 1), x = c(8, 9), dims = c(6, 1))
  sv_row <- Matrix::sparseMatrix(i = c(1, 1), j = c(3, 4), x = c(8, 9), dims = c(1, 6))

  checkEquals(x, test_StaticVector_dbl_len3(x), msg = "wrap/as double StaticVector")
  checkEquals(z, test_StaticVector_cpl_len3(z), msg = "wrap/as complex StaticVector")
  checkEquals(x, test_HybridVector_dbl_len3(x), msg = "wrap/as double HybridVector")
  checkEquals(z, test_HybridVector_cpl_len3(z), msg = "wrap/as complex HybridVector")
  checkEquals(x, test_DynamicVector_dbl(x), msg = "wrap/as double DynamicVector")
  checkEquals(z, test_DynamicVector_cpl(z), msg = "wrap/as complex DynamicVector")

  checkEquals(x, test_CustomVector1_dbl(x), msg = "wrap/as double CustomVector type 1")
  checkEquals(x, test_CustomVector2_dbl(x), msg = "wrap/as double CustomVector type 2")
  checkEquals(x, test_CustomVector3_dbl(x), msg = "wrap/as double CustomVector type 3")
  checkEquals(x, test_CustomVector4_dbl(x), msg = "wrap/as double CustomVector type 4")
  checkEquals(z, test_CustomVector1_cpl(z), msg = "wrap/as complex CustomVector type 1")
  checkEquals(z, test_CustomVector2_cpl(z), msg = "wrap/as complex CustomVector type 2")
  checkEquals(z, test_CustomVector3_cpl(z), msg = "wrap/as complex CustomVector type 3")
  checkEquals(z, test_CustomVector4_cpl(z), msg = "wrap/as complex CustomVector type 4")

  checkEquals(sv_col, test_CompressedVector_dbl_col(sv_col), msg = "wrap/as column CompressedVector")
  checkEquals(sv_row, test_CompressedVector_dbl_row(sv_row), msg = "wrap/as row CompressedVector")
}

