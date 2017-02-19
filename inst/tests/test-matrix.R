Rcpp::sourceCpp("cpp/test-matrices.cpp")
x <- matrix(c(rep(0, 4), 1, 3, rep(0, 3), 4, rep(0, 2), 2, 0, 5, 0), 4)
z <- matrix(c(rep(0, 4), 1+2i, 3, rep(0, 3), 4-1i, rep(0, 2), 2+1i, 0, 5-2i, 0), 4)

context("test - wrap/as StaticMatrix")
test_that("test - wrap/as StaticVector", {
  expect_equal(test_StaticMatrix_dbl_dim44(x), x)
  expect_equal(test_StaticMatrix_cpl_dim44(z), z)
})

context("test - wrap/as HybridMatrix")
test_that("test - wrap/as HybridMatrix", {
  expect_equal(test_HybridMatrix_dbl_dim44(x), x)
  expect_equal(test_HybridMatrix_cpl_dim44(z), z)
})

context("test - wrap/as DynamicMatrix")
test_that("test - wrap/as DynamicMatrix", {
  expect_equal(test_DynamicMatrix_dbl(x), x)
  expect_equal(test_DynamicMatrix_cpl(z), z)
})

context("test - wrap/as CustomMatrix")
test_that("test - wrap/as CustomMatrix", {
  expect_equal(test_CustomMatrix1_dbl(x), x)
  expect_equal(test_CustomMatrix2_dbl(x), x)
  expect_equal(test_CustomMatrix3_dbl(x), x)
  expect_equal(test_CustomMatrix4_dbl(x), x)
  expect_equal(test_CustomMatrix1_cpl(z), z)
  expect_equal(test_CustomMatrix2_cpl(z), z)
  expect_equal(test_CustomMatrix3_cpl(z), z)
  expect_equal(test_CustomMatrix4_cpl(z), z)
})

context("test - wrap/as CompressedMatrix")
x_sp_dgC <- as(x, "dgCMatrix")
x_sp_dgR <- as(x, "dgRMatrix")
test_that("test - wrap/as DynamicVector", {
  expect_equal(test_CompressedMatrix_dbl_dgC(x_sp_dgC), x_sp_dgC)
  expect_equal(test_CompressedMatrix_dbl_dgR(x_sp_dgC), x_sp_dgR)
})
