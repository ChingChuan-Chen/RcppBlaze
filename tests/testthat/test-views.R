Rcpp::sourceCpp("cpp/test-views.cpp")
x <- matrix(c(rep(0, 4), 1, 3, rep(0, 3), 4, rep(0, 2), 2, 0, 5, 0), 4)
x_sp <- as(x, "dgCMatrix")
z <- matrix(c(rep(0, 4), 1+2i, 3, rep(0, 3), 4-1i, rep(0, 2), 2+1i, 0, 5-2i, 0), 4)

context("test - wrap DenseColumn")
test_that("test - wrap DenseColumn", {
  expect_equal(test_column_dbl(x, 0L), x[ , 1])
  expect_equal(test_column_dbl(x, 1L), x[ , 2])
  expect_equal(test_column_dbl(x, 2L), x[ , 3])
  expect_equal(test_column_dbl(x, 3L), x[ , 4])
  expect_equal(test_column_cpl(z, 0L), z[ , 1])
  expect_equal(test_column_cpl(z, 1L), z[ , 2])
  expect_equal(test_column_cpl(z, 2L), z[ , 3])
  expect_equal(test_column_cpl(z, 3L), z[ , 4])
})

context("test - wrap DenseRow")
test_that("test - wrap DenseRow", {
  expect_equal(test_row_dbl(x, 0L), x[1, ])
  expect_equal(test_row_dbl(x, 1L), x[2, ])
  expect_equal(test_row_dbl(x, 2L), x[3, ])
  expect_equal(test_row_dbl(x, 3L), x[4, ])
  expect_equal(test_row_cpl(z, 0L), z[1, ])
  expect_equal(test_row_cpl(z, 1L), z[2, ])
  expect_equal(test_row_cpl(z, 2L), z[3, ])
  expect_equal(test_row_cpl(z, 3L), z[4, ])
})

context("test - wrap DenseSubvector")
test_that("test - wrap DenseSubvector", {
  expect_equal(test_subvec_col_dbl(x[ , 1], 0L, 2L), x[1:2, 1])
  expect_equal(test_subvec_col_dbl(x[ , 2], 0L, 2L), x[1:2, 2])
  expect_equal(test_subvec_col_dbl(x[ , 3], 0L, 2L), x[1:2, 3])
  expect_equal(test_subvec_col_dbl(x[ , 4], 0L, 2L), x[1:2, 4])
  expect_equal(test_subvec_row_dbl(x[1, ], 0L, 2L), x[1, 1:2])
  expect_equal(test_subvec_row_dbl(x[2, ], 0L, 2L), x[2, 1:2])
  expect_equal(test_subvec_row_dbl(x[3, ], 0L, 2L), x[3, 1:2])
  expect_equal(test_subvec_row_dbl(x[4, ], 0L, 2L), x[4, 1:2])

  expect_equal(test_subvec_col_cpl(z[ , 1], 0L, 2L), z[1:2, 1])
  expect_equal(test_subvec_col_cpl(z[ , 2], 0L, 2L), z[1:2, 2])
  expect_equal(test_subvec_col_cpl(z[ , 3], 0L, 2L), z[1:2, 3])
  expect_equal(test_subvec_col_cpl(z[ , 4], 0L, 2L), z[1:2, 4])
  expect_equal(test_subvec_row_cpl(z[1, ], 0L, 2L), z[1, 1:2])
  expect_equal(test_subvec_row_cpl(z[2, ], 0L, 2L), z[2, 1:2])
  expect_equal(test_subvec_row_cpl(z[3, ], 0L, 2L), z[3, 1:2])
  expect_equal(test_subvec_row_cpl(z[4, ], 0L, 2L), z[4, 1:2])
})

context("test - wrap DenseSubmatrix")
test_that("test - wrap DenseSubmatrix", {
  expect_equal(test_submat_dbl(x, 0L, 0L, 2L, 2L), x[1:2, 1:2])
  expect_equal(test_submat_cpl(z, 0L, 0L, 2L, 2L), z[1:2, 1:2])
})

context("test - wrap SparseColumn")
test_that("test - wrap SparseColumn", {
  expect_equal(test_sp_column_dbl(x_sp, 1), as(matrix(x_sp[, 2], 4), "dgCMatrix"))
  expect_equal(test_sp_column_dbl(x_sp, 2), as(matrix(x_sp[, 3], 4), "dgCMatrix"))
  expect_equal(test_sp_column_dbl(x_sp, 3), as(matrix(x_sp[, 4], 4), "dgCMatrix"))
})

context("test - wrap SparseRow")
test_that("test - wrap SparseRow", {
  expect_equal(test_sp_row_dbl(x_sp, 0), as(matrix(x_sp[1, ], 1), "dgCMatrix"))
  expect_equal(test_sp_row_dbl(x_sp, 1), as(matrix(x_sp[2, ], 1), "dgCMatrix"))
  expect_equal(test_sp_row_dbl(x_sp, 2), as(matrix(x_sp[3, ], 1), "dgCMatrix"))
})

context("test - wrap SparseSubvector")
x_colvec_sp <- as(matrix(c(0, 1, 2), 3), "dgCMatrix")
x_rowvec_sp <- as(matrix(c(0, 1, 2), 1), "dgCMatrix")
test_that("test - wrap SparseSubvector", {
  expect_equal(test_sp_subcolvec_dbl(x_colvec_sp, 0, 2), as(matrix(x_colvec_sp[1:2], 2), "dgCMatrix"))
  expect_equal(test_sp_subrowvec_dbl(x_rowvec_sp, 0, 2), as(matrix(x_rowvec_sp[1:2], 1), "dgCMatrix"))
})

context("test - wrap SparseSubmatrix")
test_that("test - wrap SparseSubmatrix", {
  expect_equal(test_sp_submat_dbl(x_sp, 0, 0, 2, 2), x_sp[1:2, 1:2])
})
