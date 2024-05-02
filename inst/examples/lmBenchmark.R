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
exprs$blaze_lu <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', X, y, 3L))

if (suppressMessages(require(RcppEigen, quietly = TRUE))) {
  # versions from RcppEigen
  # versions which can handle rank-deficient cases.
  ## column-pivoted QR decomposition - similar to lm.fit
  exprs$eigen_PivQR <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 0L))
  ## LDLt Cholesky decomposition with rank detection
  exprs$eigen_LDLt <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 2L))
  ## SVD (the JacobiSVD class from Eigen)
  exprs$eigen_SVD <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 4L))
  ## eigenvalues and eigenvectors of X'X
  exprs$eigen_SymmEig <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 5L))
  ## SVD using the Lapack subroutine dgesdd (SVD) and Eigen support
  exprs$lapack_GESDD <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 6L))

  # versions which cannot handle rank-deficient cases.
  ## Unpivoted  QR decomposition
  exprs$eigen_QR <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 1L))
  ## LLt Cholesky decomposition
  exprs$eigen_LLt <- expression(.Call("_RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", X, y, 3L))
}

if (suppressMessages(require(RcppArmadillo, quietly = TRUE))) {
  # versions from RcppArmadillo
  code <- '
  // [[Rcpp::depends(RcppArmadillo)]]
  #include <RcppArmadillo.h>
  using Rcpp::_;

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_direct(const arma::mat& X, const arma::vec& y) {
    arma::mat R = chol(X.t() * X);
    arma::vec coef = arma::solve(R, arma::solve(R.t(), X.t() * y));
    arma::vec fitted = X*coef;
    arma::vec res  = y - fitted;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::vec se = arma::sqrt(s2 * arma::sum(arma::square(arma::inv(R)), 1));
    return Rcpp::List::create(
      _["coefficients"]  = coef,
      _["stderr"]        = se,
      _["df.residual"]   = df,
      _["s"]             = std::sqrt(s2),
      _["fitted.values"] = fitted
    );
  }

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_qr(const arma::mat& X, const arma::vec& y) {
    arma::mat Q, R;
    arma::qr_econ(Q, R, X);
    arma::vec coef = arma::solve(R, Q.t() * y);
    arma::vec fitted = X*coef;
    arma::vec res  = y - fitted;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::vec se = arma::sqrt(s2 * arma::sum(arma::square(arma::inv(R)), 1));
    return Rcpp::List::create(
      _["coefficients"]  = coef,
      _["stderr"]        = se,
      _["df.residual"]   = df,
      _["s"]             = std::sqrt(s2),
      _["fitted.values"] = fitted
    );
  }

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_chol(const arma::mat& X, const arma::vec& y) {
    arma::mat xtx = X.t() * X;
    arma::vec coef = arma::solve(xtx, X.t() * y);
    arma::vec fitted = X*coef;
    arma::vec res  = y - fitted;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::colvec se = arma::sqrt(s2 * arma::diagvec(arma::inv_sympd(xtx)));
    return Rcpp::List::create(
      _["coefficients"]  = coef,
      _["stderr"]        = se,
      _["df.residual"]   = df,
      _["s"]             = std::sqrt(s2),
      _["fitted.values"] = fitted
    );
  }

  // [[Rcpp::export]]
  Rcpp::List arma_fastLm_pinv(const arma::mat& X, const arma::vec& y) {
    arma::mat xtx_inv = arma::pinv(X.t() * X);
    arma::vec coef = xtx_inv * X.t() * y;
    arma::vec fitted = X*coef;
    arma::vec res  = y - fitted;
    arma::uword df = X.n_rows - X.n_cols;
    double s2 = arma::dot(res, res) / (double) df;
    arma::colvec se = arma::sqrt(s2 * arma::diagvec(xtx_inv));
    return Rcpp::List::create(
      _["coefficients"]  = coef,
      _["stderr"]        = se,
      _["df.residual"]   = df,
      _["s"]             = std::sqrt(s2),
      _["fitted.values"] = fitted
    );
  }'
  Rcpp::sourceCpp(code = code)

  # versions  which can handle rank-deficient cases.
  ## use arma::solve to solve linear equation which uses QR decomposition
  exprs$arma_fastLm <- expression(.Call("_RcppArmadillo_fastLm_impl", PACKAGE = "RcppArmadillo", X, y))
  exprs$arma_direct_solve <- expression(arma_fastLm_direct(X, y))
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
  exprs$GSL <- expression(RcppGSL:::fastLm(X, y))
}

do_bench <- function(n = 1e4L, p = 100L, nrep = 20L, suppressSVD = (n > 1e5L || p > 2e2L)) {
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), ncol = p - 1L))
  y <- X %*% rnorm(p, sd = 3) + rnorm(n, sd = 5)
  if (suppressSVD) exprs <- exprs[!names(exprs) %in% c("eigen_SVD", "GSL")]
  cat("lm benchmark for n = ", n, " and p = ", p, ": nrep = ", nrep, "\n", sep="")
  microbenchmark(list = do.call(c, exprs), times = nrep)
}

print(do_bench())

sessionInfo()
# RcppGSL_0.3.13 RcppArmadillo_0.12.8.1.0 RcppEigen_0.3.4.0.0 microbenchmark_1.4.10 RcppBlaze_1.0.0

.Call("_RcppBlaze_blaze_version", FALSE, PACKAGE = "RcppBlaze")
