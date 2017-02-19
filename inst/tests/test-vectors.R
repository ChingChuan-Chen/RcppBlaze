Rcpp::sourceCpp("cpp/test-vectors.cpp")
x <- 1:3
z <- 1:3 + 2i

context("test - wrap/as StaticVector")
test_that("test - wrap/as StaticVector", {
  expect_equal(test_StaticVector_dbl_len3(x), x)
  expect_equal(test_StaticVector_cpl_len3(z), z)
})

context("test - wrap/as HybridVector")
test_that("test - wrap/as HybridVector", {
  expect_equal(test_HybridVector_dbl_len3(x), x)
  expect_equal(test_HybridVector_cpl_len3(z), z)
})

context("test - wrap/as DynamicVector")
test_that("test - wrap/as DynamicVector", {
  expect_equal(test_DynamicVector_dbl(x), x)
  expect_equal(test_DynamicVector_cpl(z), z)
})

context("test - wrap/as CustomVector")
test_that("test - wrap/as CustomVector", {
  expect_equal(test_CustomVector1_dbl(x), x)
  expect_equal(test_CustomVector2_dbl(x), x)
  expect_equal(test_CustomVector3_dbl(x), x)
  expect_equal(test_CustomVector4_dbl(x), x)
  expect_equal(test_CustomVector1_cpl(z), z)
  expect_equal(test_CustomVector2_cpl(z), z)
  expect_equal(test_CustomVector3_cpl(z), z)
  expect_equal(test_CustomVector4_cpl(z), z)
})

context("test - wrap/as CompressedVector")
sv_col <- Matrix::sparseMatrix(i = c(3, 4), j = c(1, 1), x = c(8, 9), dims = c(6, 1))
sv_row <- Matrix::sparseMatrix(i = c(1, 1), j = c(3, 4), x = c(8, 9), dims = c(1, 6))
test_that("test - wrap/as HybridVector", {
  expect_equal(test_CompressedVector_dbl_col(sv_col), sv_col)
  expect_equal(test_CompressedVector_dbl_row(sv_row), sv_row)
})
