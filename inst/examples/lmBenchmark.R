## lmBenchmark.R: Benchmark different implementations of linear model solutions
##
## Copyright (C)  2017 Douglas Bates, Dirk Eddelbuettel, Romain Francois and Chingchuan Chen
##
## This file is based on lmBenchmark.R from RcppEigen.
## This file is part of RcppBlaze.

require("stats", character=TRUE, quietly=TRUE)
require("microbenchmark", character=TRUE, quietly=TRUE)
require("RcppBlaze", character=TRUE, quietly=TRUE)

## define different versions of lm
exprs <- list()

# default version used in lm()
exprs$lm.fit <- expression(stats::lm.fit(mm, y))

exprs$blaze_qr <- expression(.Call('RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', mm, y, 0L))
exprs$blaze_ldlt <- expression(.Call('RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', mm, y, 1L))
exprs$blaze_llt <- expression(.Call('RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', mm, y, 2L))

if (suppressMessages(require("RcppEigen", character = TRUE, quietly = TRUE))) {
  # versions from RcppEigen
  # versions which can handle rank-deficient cases.
  ## column-pivoted QR decomposition - similar to lm.fit
  exprs$eigen_PivQR <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 0L))
  ## LDLt Cholesky decomposition with rank detection
  exprs$eigen_LDLt <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 2L))
  ## SVD (the JacobiSVD class from Eigen)
  exprs$eigen_SVD <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 4L))
  ## eigenvalues and eigenvectors of X'X
  exprs$eigen_SymmEig <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 5L))
  ## SVD using the Lapack subroutine dgesdd (SVD) and Eigen support
  exprs$lapack_GESDD <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 6L))

  # versions which cannot handle rank-deficient cases.
  ## Unpivoted  QR decomposition
  exprs$eigen_QR <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 1L))
  ## LLt Cholesky decomposition
  exprs$eigen_LLt <- expression(.Call("RcppEigen_fastLm_Impl", PACKAGE = "RcppEigen", mm, y, 3L))
}

if (suppressMessages(require("RcppArmadillo", character = TRUE, quietly = TRUE))) {
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
  exprs$arma_solve1 <- expression(.Call("RcppArmadillo_fastLm", PACKAGE = "RcppArmadillo", mm, y))
  exprs$arma_solve2 <- expression(arma_fastLm_direct(mm, y))
  ## use cholesky decomposition to solve linear equation
  exprs$arma_qr <- expression(arma_fastLm_qr(mm, y))
  ## use cholesky decomposition to solve linear equation
  exprs$arma_chol <- expression(arma_fastLm_chol(mm, y))

  # versions which can handle rank-deficient cases.
  ## use arma::solve to solve linear equation which uses LU decomposition
  exprs$arma_pinv <- expression(arma_fastLm_pinv(mm, y))
}

if (suppressMessages(require("RcppGSL", character = TRUE, quietly = TRUE))) {
  # versions from RcppGSL  (it cannot handle rank-deficient cases.)
  exprs$GSL <- expression(.Call("RcppGSL_fastLm", PACKAGE = "RcppGSL", mm, y))
}

do_bench <- function(n = 1e5L, p = 40L, nrep = 20L, suppressSVD = (n > 1e5L || p > 2e2L)) {
  mm <- cbind(1, matrix(rnorm(n * (p - 1L)), nc = p-1L))
  y <- mm %*% rnorm(p, sd = 3) + rnorm(n, sd = 5)
  if (suppressSVD) exprs <- exprs[!names(exprs) %in% c("eigen_SVD", "GSL")]
  cat("lm benchmark for n = ", n, " and p = ", p, ": nrep = ", nrep, "\n", sep="")
  microbenchmark(list = do.call(c, exprs), times = nrep)
}

print(do_bench())
# Reference benchmark:
## Unit: milliseconds
##           expr         min          lq        mean      median          uq        max neval
##         lm.fit  115.051560  120.155425  128.755254  127.093843  132.403677  183.90210    20
##       blaze_qr  106.871918  112.214372  116.391643  113.883625  119.122072  133.05420    20
##     blaze_ldlt   42.476913   43.073899   45.385867   44.110764   47.243592   53.39617    20
##      blaze_llt   42.067316   42.787768   48.667232   44.640606   48.425277  106.40351    20
##    eigen_PivQR  116.095153  117.692579  124.490097  121.240856  125.342379  184.68326    20
##     eigen_LDLt   20.843782   21.379182   23.853525   22.242114   27.287612   31.50733    20
##      eigen_SVD  905.852198  913.107176  925.425030  916.658670  935.456956  958.27382    20
##  eigen_SymmEig   54.920454   56.233064   62.434809   61.254718   66.790999   77.00092    20
##   lapack_GESDD  135.296744  141.577468  150.809217  152.081570  158.851469  167.68764    20
##       eigen_QR  111.750797  113.755041  116.631769  117.010017  118.990123  123.88158    20
##      eigen_LLt   20.859288   21.564671   23.103130   21.729387   22.806041   30.04537    20
##    arma_solve1   57.162995   58.569227   63.041392   61.149393   66.653785   75.89618    20
##    arma_solve2    7.276483    9.572271    9.411490    9.730258    9.821393   10.15975    20
##        arma_qr   81.104203   85.088260   88.977305   89.132587   90.633319  105.91141    20
##      arma_chol    7.115569    8.129028    9.073557    9.607965    9.819344   10.31159    20
##      arma_pinv    7.434177    9.764927   10.151937   10.124933   10.418380   12.56379    20
##            GSL 2331.238345 2380.913915 2475.658602 2469.485940 2558.890029 2628.49481    20

print(do_bench(n = 1.6e4L, p = 250L))
# Reference benchmark:
## Unit: milliseconds
##           expr       min        lq      mean    median        uq       max neval
##         lm.fit 407.83260 419.55832 432.62875 427.20578 440.99923 478.96781    20
##       blaze_qr 219.92874 227.78510 239.28999 237.47001 250.19207 266.19822    20
##     blaze_ldlt  58.72502  60.80884  70.47990  70.92602  74.01672  91.23937    20
##      blaze_llt  57.32625  61.57011  67.55130  66.95323  73.42149  80.17939    20
##    eigen_PivQR 518.99621 522.27196 526.01878 524.51128 526.43902 543.24842    20
##     eigen_LDLt  75.36693  76.20163  88.25299  81.12176  93.58402 156.64169    20
##  eigen_SymmEig 238.34362 239.23668 252.46615 246.19090 264.70773 286.33691    20
##   lapack_GESDD 379.06839 390.20912 403.39123 399.75418 412.68441 437.22115    20
##       eigen_QR 229.67041 236.57431 245.55035 242.44104 252.59377 277.47120    20
##      eigen_LLt  77.88273  78.69592  84.84525  81.22255  85.03238 106.71774    20
##    arma_solve1 109.92370 111.43980 120.32004 118.31751 129.97521 136.11623    20
##    arma_solve2  18.69223  19.28234  24.04288  22.62406  29.70716  31.36485    20
##        arma_qr 200.53054 209.84492 220.59700 213.80571 231.85517 252.98640    20
##      arma_chol  17.30194  17.94735  21.63480  18.82608  28.22968  29.13226    20
##      arma_pinv  28.12055  31.52635  36.17278  38.00544  40.82814  42.26187    20

print(do_bench(n = 50, p = 10))
# Reference benchmark:
## Unit: microseconds
##           expr    min       lq      mean   median       uq     max neval
##         lm.fit 71.095 102.1075 112.18640 107.8130 124.7810 191.341    20
##       blaze_qr 19.310  26.6245  29.79880  30.5740  33.5000  38.327    20
##     blaze_ldlt 11.704  15.5060  18.81285  18.5790  20.9190  31.891    20
##      blaze_llt 13.167  14.9220  18.18365  17.9930  21.0655  26.624    20
##    eigen_PivQR 39.205  44.9105  57.82705  56.4665  69.1935  92.746    20
##     eigen_LDLt 36.280  37.8885  54.08205  54.7110  65.6820  80.458    20
##      eigen_SVD 76.947  81.7740 106.59845 117.1750 124.1970 133.412    20
##  eigen_SymmEig 46.519  55.4425  67.35000  71.0945  76.6540  85.138    20
##   lapack_GESDD 84.553 119.5150 121.53400 122.7335 134.2905 145.408    20
##       eigen_QR 41.546  44.1785  61.49865  57.6365  67.1455 162.376    20
##      eigen_LLt 39.790  44.6175  55.12070  56.6130  63.7810  78.409    20
##    arma_solve1 92.745  95.8175 131.48125 141.7505 157.8415 198.070    20
##    arma_solve2 30.135  34.0855  41.45770  42.1305  48.7135  59.100    20
##        arma_qr 25.454  32.6220  45.59745  41.1065  45.7875 105.325    20
##      arma_chol 24.576  30.5740  35.21140  36.1330  38.3270  46.811    20
##      arma_pinv 44.178  47.6900  60.97205  65.3895  70.0710  88.941    20
##            GSL 42.130  45.0565  56.96370  60.4155  64.8040  71.095    20

sessionInfo()
# Reference benchmark runs on windows 10 x64 with MRO-3.3.1,
# RcppGSL 0.3.0, RcppArmadillo 0.7.600.1.0 and RcppEigen 0.3.2.9.0
# and CPU is i7-3770K@4.2GHz.

.Call("RcppBlaze_blaze_version", FALSE, PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_SSE", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_AVX", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_AVX2", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_MIC", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_FMA", PACKAGE = "RcppBlaze")
