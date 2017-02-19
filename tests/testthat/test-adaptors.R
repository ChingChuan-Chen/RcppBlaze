Sys.setenv(R_TEST = "")
Rcpp::sourceCpp("cpp/test-adaptors.cpp")
x <- matrix(c(rep(0, 4), 1, 3, rep(0, 3), 4, rep(0, 2), 2, 0, 5, 0), 4)
x2 <- t(x)
z <- matrix(c(rep(0, 4), 1+2i, 3, rep(0, 3), 4-1i, rep(0, 2), 2+1i, 0, 5-2i, 0), 4)
z2 <- t(z)

context("test - wrap DiagonalMatrix")
test_that("test - wrap DiagonalMatrix", {
  expect_equal(test_diag_matrix_dbl(x), diag(diag(x)))
  expect_equal(test_diag_matrix_cpl(z), diag(diag(z)))
})

context("test - wrap LowerMatrix")
test_that("test - wrap LowerMatrix", {
  expect_equal(test_lower_matrix_dbl(x2), x2)
  expect_equal(test_lower_matrix_cpl(z2), z2)
})

context("test - wrap UpperMatrix")
test_that("test - wrap UpperMatrix", {
  expect_equal(test_upper_matrix_dbl(x), x)
  expect_equal(test_upper_matrix_cpl(z), z)
})

context("test - wrap StrictlyLowerMatrix")
test_that("test - wrap StrictlyLowerMatrix", {
  expect_equal(test_strictlylower_matrix_dbl(x2), x2 - diag(diag(x2)))
  expect_equal(test_strictlylower_matrix_cpl(z2), z2 - diag(diag(z2)))
})

context("test - wrap StrictlyUpperMatrix")
test_that("test - wrap StrictlyUpperMatrix",{
  expect_equal(test_strictlyupper_matrix_dbl(x), x - diag(diag(x)))
  expect_equal(test_strictlyupper_matrix_cpl(z), z - diag(diag(z)))
})

context("test - wrap UniLowerMatrix")
test_that("test - wrap UniLowerMatrix", {
  expect_equal(test_unilower_matrix_dbl(x2), x2 - diag(diag(x2)) + diag(1, nrow(x2)))
  expect_equal(test_unilower_matrix_cpl(z2), z2 - diag(diag(z2)) + diag(1, nrow(z2)))
})

context("test - wrap UniUpperMatrix")
test_that("test - wrap UniUpperMatrix",{
  expect_equal(test_uniupper_matrix_dbl(x), x - diag(diag(x)) + diag(1, nrow(x)))
  expect_equal(test_uniupper_matrix_cpl(z), z - diag(diag(z)) + diag(1, nrow(z)))
})

context("test - wrap SymmetricMatrix")
test_that("test - wrap SymmetricMatrix", {
  expect_equal(test_symmetric_matrix(x), x + t(x) - diag(diag(x)))
})

context("test - wrap HermitianMatrix")
test_that("test - wrap HermitianMatrix", {
  expect_equal(test_hermitian_matrix(z), z + t(Conj(z)) - diag(diag(z)))
})
