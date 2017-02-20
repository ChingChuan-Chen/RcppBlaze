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

cppFile <- "test-views.cpp"
if (file.exists(file.path("cpp", cppFile))) {
  sourceCpp(file.path( "cpp", cppFile))
} else {
  sourceCpp(system.file("unitTests", "cpp", cppFile, package = "RcppBlaze"))
}

testViews <- function(){
  x <- matrix(c(rep(0, 4), 1, 3, rep(0, 3), 4, rep(0, 2), 2, 0, 5, 0), 4)
  x_sp <- as(x, "dgCMatrix")
  z <- matrix(c(rep(0, 4), 1+2i, 3, rep(0, 3), 4-1i, rep(0, 2), 2+1i, 0, 5-2i, 0), 4)

  checkEquals(x[ , 1], test_column_dbl(x, 0L), msg = "wrap/as double DenseColumn 1")
  checkEquals(x[ , 2], test_column_dbl(x, 1L), msg = "wrap/as double DenseColumn 2")
  checkEquals(x[ , 3], test_column_dbl(x, 2L), msg = "wrap/as double DenseColumn 3")
  checkEquals(x[ , 4], test_column_dbl(x, 3L), msg = "wrap/as double DenseColumn 4")

  checkEquals(z[ , 1], test_column_cpl(z, 0L), msg = "wrap/as complex DenseColumn 1")
  checkEquals(z[ , 2], test_column_cpl(z, 1L), msg = "wrap/as complex DenseColumn 2")
  checkEquals(z[ , 3], test_column_cpl(z, 2L), msg = "wrap/as complex DenseColumn 3")
  checkEquals(z[ , 4], test_column_cpl(z, 3L), msg = "wrap/as complex DenseColumn 4")

  checkEquals(x[1, ], test_row_dbl(x, 0L), msg = "wrap/as double DenseRow 1")
  checkEquals(x[2, ], test_row_dbl(x, 1L), msg = "wrap/as double DenseRow 2")
  checkEquals(x[3, ], test_row_dbl(x, 2L), msg = "wrap/as double DenseRow 3")
  checkEquals(x[4, ], test_row_dbl(x, 3L), msg = "wrap/as double DenseRow 4")

  checkEquals(z[1, ], test_row_cpl(z, 0L), msg = "wrap/as complex DenseRow 1")
  checkEquals(z[2, ], test_row_cpl(z, 1L), msg = "wrap/as complex DenseRow 2")
  checkEquals(z[3, ], test_row_cpl(z, 2L), msg = "wrap/as complex DenseRow 3")
  checkEquals(z[4, ], test_row_cpl(z, 3L), msg = "wrap/as complex DenseRow 4")

  checkEquals(x[1:2, 1], test_subvec_col_dbl(x[ , 1], 0L, 2L), msg = "wrap/as double Column DenseSubvector 1")
  checkEquals(x[1:2, 2], test_subvec_col_dbl(x[ , 2], 0L, 2L), msg = "wrap/as double Column DenseSubvector 2")
  checkEquals(x[1:2, 3], test_subvec_col_dbl(x[ , 3], 0L, 2L), msg = "wrap/as double Column DenseSubvector 3")
  checkEquals(x[1:2, 4], test_subvec_col_dbl(x[ , 4], 0L, 2L), msg = "wrap/as double Column DenseSubvector 4")
  checkEquals(z[1:2, 1], test_subvec_col_cpl(z[ , 1], 0L, 2L), msg = "wrap/as complex Column DenseSubvector 1")
  checkEquals(z[1:2, 2], test_subvec_col_cpl(z[ , 2], 0L, 2L), msg = "wrap/as complex Column DenseSubvector 2")
  checkEquals(z[1:2, 3], test_subvec_col_cpl(z[ , 3], 0L, 2L), msg = "wrap/as complex Column DenseSubvector 3")
  checkEquals(z[1:2, 4], test_subvec_col_cpl(z[ , 4], 0L, 2L), msg = "wrap/as complex Column DenseSubvector 4")

  checkEquals(x[1, 1:2], test_subvec_col_dbl(x[1 , ], 0L, 2L), msg = "wrap/as double Row DenseSubvector 1")
  checkEquals(x[2, 1:2], test_subvec_col_dbl(x[2 , ], 0L, 2L), msg = "wrap/as double Row DenseSubvector 2")
  checkEquals(x[3, 1:2], test_subvec_col_dbl(x[3 , ], 0L, 2L), msg = "wrap/as double Row DenseSubvector 3")
  checkEquals(x[4, 1:2], test_subvec_col_dbl(x[4 , ], 0L, 2L), msg = "wrap/as double Row DenseSubvector 4")
  checkEquals(z[1, 1:2], test_subvec_col_cpl(z[1 , ], 0L, 2L), msg = "wrap/as complex Row DenseSubvector 1")
  checkEquals(z[2, 1:2], test_subvec_col_cpl(z[2 , ], 0L, 2L), msg = "wrap/as complex Row DenseSubvector 2")
  checkEquals(z[3, 1:2], test_subvec_col_cpl(z[3 , ], 0L, 2L), msg = "wrap/as complex Row DenseSubvector 3")
  checkEquals(z[4, 1:2], test_subvec_col_cpl(z[4 , ], 0L, 2L), msg = "wrap/as complex Row DenseSubvector 4")

  checkEquals(x[1:2, 1:2], test_submat_dbl(x, 0L, 0L, 2L, 2L), msg = "wrap/as double DenseSubmatrix")
  checkEquals(z[1:2, 1:2], test_submat_cpl(z, 0L, 0L, 2L, 2L), msg = "wrap/as complex DenseSubmatrix")

  checkEquals(as(matrix(x_sp[, 2], 4), "dgCMatrix"), test_sp_column_dbl(x_sp, 1), msg = "wrap/as double SparseColumn 1")
  checkEquals(as(matrix(x_sp[, 3], 4), "dgCMatrix"), test_sp_column_dbl(x_sp, 2), msg = "wrap/as double SparseColumn 2")
  checkEquals(as(matrix(x_sp[, 4], 4), "dgCMatrix"), test_sp_column_dbl(x_sp, 3), msg = "wrap/as double SparseColumn 3")

  checkEquals(as(matrix(x_sp[1 , ], 1), "dgCMatrix"), test_sp_row_dbl(x_sp, 0), msg = "wrap/as double SparseRow 1")
  checkEquals(as(matrix(x_sp[2 , ], 1), "dgCMatrix"), test_sp_row_dbl(x_sp, 1), msg = "wrap/as double SparseRow 2")
  checkEquals(as(matrix(x_sp[3 , ], 1), "dgCMatrix"), test_sp_row_dbl(x_sp, 2), msg = "wrap/as double SparseRow 3")

  x_colvec_sp <- as(matrix(c(0, 1, 2), 3), "dgCMatrix")
  x_rowvec_sp <- as(matrix(c(0, 1, 2), 1), "dgCMatrix")
  checkEquals(as(matrix(x_colvec_sp[1:2], 2), "dgCMatrix"), test_sp_subcolvec_dbl(x_colvec_sp, 0, 2),
              msg = "wrap/as double Column SparseSubvector")
  checkEquals(as(matrix(x_rowvec_sp[1:2], 1), "dgCMatrix"), test_sp_subrowvec_dbl(x_rowvec_sp, 0, 2),
              msg = "wrap/as double Row SparseSubvector")

  checkEquals(x_sp[1:2, 1:2], test_sp_submat_dbl(x_sp, 0, 0, 2, 2), msg = "wrap/as double SparseSubmatrix")
}
