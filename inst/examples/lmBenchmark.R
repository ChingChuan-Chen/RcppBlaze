## lmBenchmark.R: Benchmark different implementations of linear model solutions
##
## Copyright (C) 2017 - 2024 Douglas Bates, Dirk Eddelbuettel, Romain Francois and Ching-Chuan Chen
##
## This file is based on lmBenchmark.R from RcppEigen.
## This file is part of RcppBlaze.

suppressPackageStartupMessages({
  require(stats)
  require(microbenchmark)
  require(RcppBlaze)
})

## define different versions of lm
exprs <- list()

# default version used in lm()
exprs$lm.fit <- expression(stats::lm.fit(X, y))

exprs$blaze_qr <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', X, y, 0L))
exprs$blaze_ldlt <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', X, y, 1L))
exprs$blaze_llt <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', X, y, 2L))

if (suppressMessages(require(RcppEigen, quietly = TRUE))) {
  # versions from RcppEigen
  # versions which can handle rank-deficient cases.
  ## column-pivoted QR decomposition - similar to lm.fit
  exprs$eigen_PivQR <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 0L))
  ## LDLt Cholesky decomposition with rank detection
  exprs$eigen_LDLt <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 2L))
  ## SVD (the JacobiSVD class from Eigen)
  exprs$eigen_SVD <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 4L))
  ## eigenvalues and eigenvectors of X'X
  exprs$eigen_SymmEig <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 5L))
  ## SVD using the Lapack subroutine dgesdd (SVD) and Eigen support
  exprs$lapack_GESDD <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 6L))

  # versions which cannot handle rank-deficient cases.
  ## Unpivoted  QR decomposition
  exprs$eigen_QR <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 1L))
  ## LLt Cholesky decomposition
  exprs$eigen_LLt <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 3L))
}

if (suppressMessages(require(RcppArmadillo, quietly = TRUE))) {
  # versions from RcppArmadillo
  code <- '
  // [[Rcpp::depends(RcppArmadillo)]]
  #include <RcppArmadillo.h>
  using Rcpp::_;
  using Rcpp::List;

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_direct(const arma::mat& X, const arma::vec& y) {
    arma::mat R = chol(X.t() * X);
    arma::vec coef = arma::solve(R, arma::solve(R.t(), X.t() * y));
    arma::vec res  = y - X*coef;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::vec se = arma::sqrt(s2 * arma::sum(arma::square(arma::inv(R)), 1));
    return List::create(_["coefficients"] = coef,
                        _["stderr"]       = se,
                        _["df.residual"]  = df);
  }

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_qr(const arma::mat& X, const arma::vec& y) {
    arma::mat Q, R;
    arma::qr_econ(Q, R, X);
    arma::vec coef = arma::solve(R, Q.t() * y);
    arma::vec res  = y - X*coef;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::vec se = arma::sqrt(s2 * arma::sum(arma::square(arma::inv(R)), 1));
    return List::create(_["coefficients"] = coef,
                        _["stderr"]       = se,
                        _["df.residual"]  = df);
  }

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_chol(const arma::mat& X, const arma::vec& y) {
    arma::mat xtx = X.t() * X;
    arma::vec coef = arma::solve(xtx, X.t() * y);
    arma::vec res  = y - X*coef;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::colvec se = arma::sqrt(s2 * arma::diagvec(arma::inv_sympd(xtx)));
    return List::create(_["coefficients"] = coef,
                        _["stderr"]       = se,
                        _["df.residual"]  = df);
  }

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_pinv(const arma::mat& X, const arma::vec& y) {
    arma::mat xtx_inv = arma::pinv(X.t() * X);
    arma::vec coef = xtx_inv * X.t() * y;
    arma::vec res  = y - X*coef;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::colvec se = arma::sqrt(s2 * arma::diagvec(xtx_inv));
    return List::create(_["coefficients"] = coef,
                        _["stderr"]       = se,
                        _["df.residual"]  = df);
  }'
  Rcpp::sourceCpp(code = code)

  # versions  which can handle rank-deficient cases.
  ## use arma::solve to solve linear equation which uses QR decomposition
  exprs$arma_solve1 <- expression(.Call("_RcppArmadillo_fastLm_impl", PACKAGE = "RcppArmadillo", X, y))
  exprs$arma_solve2 <- expression(arma_fastLm_direct(X, y))
  ## use cholesky decomposition to solve linear equation
  exprs$arma_qr <- expression(arma_fastLm_qr(X, y))
  ## use cholesky decomposition to solve linear equation
  exprs$arma_chol <- expression(arma_fastLm_chol(X, y))

  # versions which can handle rank-deficient cases.
  ## use arma::solve to solve linear equation which uses LU decomposition
  exprs$arma_pinv <- expression(arma_fastLm_pinv(X, y))
}

if (suppressMessages(require(RcppGSL, quietly = TRUE))) {
  # versions from RcppGSL  (it cannot handle rank-deficient cases.)
  exprs$GSL <- expression(.Call("RcppGSL_fastLm", PACKAGE = "RcppGSL", X, y))
}

do_bench <- function(n = 3e2L, p = 30L, nrep = 20L, suppressSVD = (n > 1e5L || p > 2e2L)) {
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), ncol = p - 1L))
  y <- X %*% rnorm(p, sd = 3) + rnorm(n, sd = 5)
  if (suppressSVD) exprs <- exprs[!names(exprs) %in% c("eigen_SVD", "GSL")]
  cat("lm benchmark for n = ", n, " and p = ", p, ": nrep = ", nrep, "\n", sep="")
  microbenchmark(list = do.call(c, exprs), times = nrep)
}

print(do_bench())

sessionInfo()
# RcppGSL_0.3.2 RcppArmadillo_0.7.960.1.2 RcppEigen_0.3.3.3.0 RcppBlaze_0.2.2

.Call("_RcppBlaze_blaze_version", FALSE, PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_SSE", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_AVX", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_AVX2", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_MIC", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_FMA", PACKAGE = "RcppBlaze")
