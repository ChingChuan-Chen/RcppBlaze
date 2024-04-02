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

# cppFile <- "test-adaptors.cpp"
# if (file.exists(file.path("cpp", cppFile))) {
#   sourceCpp(file.path( "cpp", cppFile))
# } else {
#   sourceCpp(system.file("unitTests", "cpp", cppFile, package = "RcppBlaze"))
# }
#
# testAdaptors <- function(){
#   x <- matrix(c(rep(0, 4), 1, 3, rep(0, 3), 4, rep(0, 2), 2, 0, 5, 0), 4)
#   x2 <- t(x)
#   z <- matrix(c(rep(0, 4), 1+2i, 3, rep(0, 3), 4-1i, rep(0, 2), 2+1i, 0, 5-2i, 0), 4)
#   z2 <- t(z)
#
#   checkEquals(diag(diag(x)), test_diag_matrix_dbl(x), msg = "wrap/as double DiagonalMatrix")
#   checkEquals(diag(diag(z)), test_diag_matrix_cpl(z), msg = "wrap/as complex DiagonalMatrix")
#
#   checkEquals(x, test_upper_matrix_dbl(x), msg = "wrap/as double UpperMatrix")
#   checkEquals(z, test_upper_matrix_cpl(z), msg = "wrap/as complex UpperMatrix")
#   checkEquals(x2, test_lower_matrix_dbl(x2), msg = "wrap/as double LowerMatrix")
#   checkEquals(z2, test_lower_matrix_cpl(z2), msg = "wrap/as complex LowerMatrix")
#
#   checkEquals(x - diag(diag(x)), test_strictlyupper_matrix_dbl(x), msg = "wrap/as double StrictlyUpperMatrix")
#   checkEquals(z - diag(diag(z)), test_strictlyupper_matrix_cpl(z), msg = "wrap/as complex StrictlyUpperMatrix")
#   checkEquals(x2 - diag(diag(x2)), test_strictlylower_matrix_dbl(x2), msg = "wrap/as double StrictlyLowerMatrix")
#   checkEquals(z2 - diag(diag(z2)), test_strictlylower_matrix_cpl(z2), msg = "wrap/as complex StrictlyLowerMatrix")
#
#   checkEquals(x - diag(diag(x)) + diag(1, nrow(x)), test_uniupper_matrix_dbl(x), msg = "wrap/as double UniUpperMatrix")
#   checkEquals(z - diag(diag(z)) + diag(1, nrow(z)), test_uniupper_matrix_cpl(z), msg = "wrap/as complex UniUpperMatrix")
#   checkEquals(x2 - diag(diag(x2)) + diag(1, nrow(x2)), test_unilower_matrix_dbl(x2), msg = "wrap/as double UniLowerMatrix")
#   checkEquals(z2 - diag(diag(z2)) + diag(1, nrow(z2)), test_unilower_matrix_cpl(z2), msg = "wrap/as complex UniLowerMatrix")
#
#   checkEquals(x + x2 - diag(diag(x)), test_symmetric_matrix(x), msg = "wrap/as SymmetricMatrix")
#   checkEquals(z + Conj(z2) - diag(diag(z)), test_hermitian_matrix(z), msg = "wrap/as HermitianMatrix")
# }
