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

# cppFile <- "test-matrices.cpp"
# if (file.exists(file.path("cpp", cppFile))) {
#   sourceCpp(file.path( "cpp", cppFile))
# } else {
#   sourceCpp(system.file("unitTests", "cpp", cppFile, package = "RcppBlaze"))
# }
#
# testMatrices <- function(){
#   x <- matrix(c(rep(0, 4), 1, 3, rep(0, 3), 4, rep(0, 2), 2, 0, 5, 0), 4)
#   z <- matrix(c(rep(0, 4), 1+2i, 3, rep(0, 3), 4-1i, rep(0, 2), 2+1i, 0, 5-2i, 0), 4)
#   x_sp_dgC <- as(x, "dgCMatrix")
#   x_sp_dgR <- as(x, "dgRMatrix")
#
#   checkEquals(x, test_StaticMatrix_dbl_dim44(x), msg = "wrap/as double StaticMatrix")
#   checkEquals(z, test_StaticMatrix_cpl_dim44(z), msg = "wrap/as complex StaticMatrix")
#   checkEquals(x, test_HybridMatrix_dbl_dim44(x), msg = "wrap/as double HybridMatrix")
#   checkEquals(z, test_HybridMatrix_cpl_dim44(z), msg = "wrap/as complex HybridMatrix")
#   checkEquals(x, test_DynamicMatrix_dbl(x), msg = "wrap/as double DynamicMatrix")
#   checkEquals(z, test_DynamicMatrix_cpl(z), msg = "wrap/as complex DynamicMatrix")
#
#   checkEquals(x, test_CustomMatrix1_dbl(x), msg = "wrap/as double CustomMatrix type 1")
#   checkEquals(x, test_CustomMatrix2_dbl(x), msg = "wrap/as double CustomMatrix type 2")
#   checkEquals(x, test_CustomMatrix3_dbl(x), msg = "wrap/as double CustomMatrix type 3")
#   checkEquals(x, test_CustomMatrix4_dbl(x), msg = "wrap/as double CustomMatrix type 4")
#   checkEquals(z, test_CustomMatrix1_cpl(z), msg = "wrap/as complex CustomMatrix type 1")
#   checkEquals(z, test_CustomMatrix2_cpl(z), msg = "wrap/as complex CustomMatrix type 2")
#   checkEquals(z, test_CustomMatrix3_cpl(z), msg = "wrap/as complex CustomMatrix type 3")
#   checkEquals(z, test_CustomMatrix4_cpl(z), msg = "wrap/as complex CustomMatrix type 4")
#
#   checkEquals(x_sp_dgC, test_CompressedMatrix_dbl_dgC(x_sp_dgC), msg = "wrap/as column-major CompressedMatrix")
#   checkEquals(x_sp_dgR, test_CompressedMatrix_dbl_dgR(x_sp_dgC), msg = "wrap/as row-major CompressedMatrix")
# }
