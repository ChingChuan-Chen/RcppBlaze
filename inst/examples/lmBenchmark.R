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

if (suppressMessages(require("RcppLAPACKE", character = TRUE, quietly = TRUE))) {
  Rcpp::sourceCpp(system.file("examples","RcppBlaze_with_RcppLAPACKE.cpp", package = "RcppBlaze"))
  # Let RcppBlaze use CBLAS
  exprs$blaze_qr_cblas <- expression(fastLmPure2(mm, y, 0L))
  ## use cholesky decomposition to solve linear equation
  exprs$blaze_ldlt_cblas <- expression(fastLmPure2(mm, y, 1L))
  ## use cholesky decomposition to solve linear equation
  exprs$arma_chol_cblas <- expression(fastLmPure2(mm, y, 2L))
}

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

  if (suppressMessages(require("RcppLAPACKE", character = TRUE, quietly = TRUE)) && packageVersion("RcppEigen") >= "0.3.3.3.0") {
    Rcpp::sourceCpp(system.file("examples","RcppEigen_with_RcppLAPACKE.cpp", package = "RcppBlaze"))
    # Let RcppEigen use LAPACKE
    exprs$eigen_PivQR_lapacke <- expression(fastLm_Impl2(mm, y, 0L))
    exprs$eigen_LDLt_lapacke <- expression(fastLm_Impl2(mm, y, 2L))
    exprs$eigen_SVD_lapacke <- expression(fastLm_Impl2(mm, y, 4L))
    exprs$eigen_SymmEig_lapacke <- expression(fastLm_Impl2(mm, y, 5L))
    exprs$eigen_QR_lapacke <- expression(fastLm_Impl2(mm, y, 1L))
    exprs$eigen_LLt_lapacke <- expression(fastLm_Impl2(mm, y, 3L))
  }
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
##                   expr         min         lq       mean     median         uq        max neval
##                 lm.fit  114.214261  119.03872  126.74593  120.95710  126.13894  215.48059    20
##               blaze_qr  118.354698  123.38308  126.63665  126.33408  130.31038  136.23242    20
##             blaze_ldlt   47.274765   50.70631   52.33932   51.84835   53.54847   60.34119    20
##              blaze_llt   46.194599   50.37936   52.31052   52.52404   54.09923   56.72123    20
##         blaze_qr_cblas  120.548381  121.35236  124.50694  124.34519  126.81316  131.67039    20
##       blaze_ldlt_cblas   48.511454   50.36049   52.67363   52.68027   54.58021   57.07963    20
##        arma_chol_cblas   48.084303   51.24273   52.14098   52.00092   53.86986   54.80257    20
##            eigen_PivQR  141.767241  143.99267  146.76714  145.24501  149.95859  154.18299    20
##             eigen_LDLt   22.342910   24.02913   25.38614   24.88826   27.02899   28.75280    20
##              eigen_SVD 1023.948822 1035.08166 1045.49101 1043.55314 1052.11312 1098.42519    20
##          eigen_SymmEig   57.388582   63.87469   67.23607   68.70486   70.77405   74.86256    20
##           lapack_GESDD  153.768128  160.42963  169.42494  161.53876  165.46752  275.36976    20
##               eigen_QR  136.474667  139.97686  142.19821  141.40343  144.22336  149.80236    20
##              eigen_LLt   21.639867   24.34452   25.03931   24.99256   25.79757   27.28030    20
##    eigen_PivQR_lapacke  115.408235  116.86845  120.13520  119.29238  122.10924  133.51826    20
##     eigen_LDLt_lapacke    9.651560   10.64688   11.37988   11.12230   11.83208   16.85198    20
##      eigen_SVD_lapacke  141.644654  145.10984  148.87550  148.67494  152.10897  157.21401    20
##  eigen_SymmEig_lapacke   26.196338   27.87876   29.96439   29.01451   30.61442   37.23701    20
##       eigen_QR_lapacke   47.335618   49.50546   51.02361   50.90979   52.87161   54.57553    20
##      eigen_LLt_lapacke    9.320956   10.46680   16.52613   11.16970   12.06365  110.35966    20
##            arma_solve1   58.857571   61.33139   63.02119   63.35553   64.07452   68.36928    20
##            arma_solve2    9.387370   10.26786   10.79375   10.60826   11.32286   13.06701    20
##                arma_qr   91.332731   93.97785   97.40196   96.84298   99.13613  109.79032    20
##              arma_chol    9.697493   10.13123   10.81310   10.63591   11.27883   14.77444    20
##              arma_pinv    9.676429   10.52298   11.07691   10.96359   11.63050   12.63167    20
##                    GSL 3480.072394 3553.28384 3606.27634 3582.86461 3631.09388 3866.51978    20

print(do_bench(n = 1.6e4L, p = 250L))
# Reference benchmark:
## Unit: milliseconds
##                   expr       min        lq      mean    median        uq       max neval
##                 lm.fit 507.19208 515.84071 531.92997 523.58955 539.57025 636.88851    20
##               blaze_qr 276.42330 281.71704 289.55407 284.46704 288.69905 360.58341    20
##             blaze_ldlt  76.18673  80.88334  85.62775  86.67006  90.48472  91.58288    20
##              blaze_llt  72.16039  79.99232  85.37211  83.62149  93.35277 100.46293    20
##         blaze_qr_cblas 275.21587 283.52044 291.20242 287.94086 292.80965 366.89792    20
##       blaze_ldlt_cblas 216.33986 220.20982 225.73885 223.68057 231.06764 245.94317    20
##        arma_chol_cblas 200.45044 216.21420 219.61089 220.90979 225.76687 232.61665    20
##            eigen_PivQR 655.79314 665.74868 681.73404 673.26317 682.55487 736.84120    20
##             eigen_LDLt  76.64430  84.77845  90.71813  92.53651  96.73955 102.88862    20
##          eigen_SymmEig 242.84106 282.08202 288.03181 291.08071 298.70053 309.55762    20
##           lapack_GESDD 471.23622 481.95332 494.20930 489.43723 510.17599 521.48407    20
##               eigen_QR 236.31179 258.78826 270.46830 276.03974 281.08904 293.70681    20
##              eigen_LLt  79.47813  92.88057  96.68897  96.70254  98.85571 140.21458    20
##    eigen_PivQR_lapacke 283.03740 288.18501 294.40697 291.30277 295.78376 328.70158    20
##     eigen_LDLt_lapacke  20.32653  23.00616  24.50646  24.22252  25.74329  31.25632    20
##      eigen_SVD_lapacke 361.96550 366.19298 371.34003 369.84600 375.64252 386.90701    20
##  eigen_SymmEig_lapacke  69.42399  71.98529  75.13241  74.01806  76.90557  84.47082    20
##       eigen_QR_lapacke  91.30260  94.56211  97.50923  96.05567 100.16612 108.25522    20
##      eigen_LLt_lapacke  22.16883  24.07902  26.48001  25.01684  26.66079  53.24259    20
##            arma_solve1 133.54634 137.99968 140.21952 138.90123 141.61496 156.91032    20
##            arma_solve2  23.66327  25.57214  26.94894  26.80678  28.33458  30.94824    20
##                arma_qr 246.62017 257.32937 266.93341 262.10614 266.82835 341.33851    20
##              arma_chol  22.16035  23.69663  24.75318  24.64016  26.30342  27.18639    20
##              arma_pinv  33.92016  35.88652  36.97149  36.26261  37.96727  41.80636    20

print(do_bench(n = 50, p = 10))
# Reference benchmark:
## Unit: microseconds
##                   expr    min      lq      mean   median       uq      max neval
##                 lm.fit 76.361 78.7020  98.84495  84.6990 105.4715  155.063    20
##               blaze_qr 17.262 21.5045  24.57640  23.4060  25.8935   35.694    20
##             blaze_ldlt 10.826 12.7275  15.97495  14.4830  18.2865   27.795    20
##              blaze_llt  9.948 14.1900  16.41370  15.6530  19.1635   26.039    20
##         blaze_qr_cblas 17.847 22.2360  26.17080  24.7225  31.1590   37.449    20
##       blaze_ldlt_cblas 11.411 15.2140  18.63715  17.7010  20.6265   33.061    20
##        arma_chol_cblas 11.411 14.7750  81.14450  16.2380  22.2360 1243.711    20
##            eigen_PivQR 42.423 45.2025  50.80510  47.1050  54.4180   71.680    20
##             eigen_LDLt 39.205 40.9600  55.93970  45.9340  60.4160  177.883    20
##              eigen_SVD 91.575 94.9390 112.58110  99.4735 136.0450  151.551    20
##          eigen_SymmEig 51.200 53.6870  60.98640  56.7585  70.8020   77.824    20
##           lapack_GESDD 82.797 85.4305  98.44995  87.0395 115.2730  136.923    20
##               eigen_QR 42.423 44.3250  51.44895  48.5670  60.8550   67.877    20
##              eigen_LLt 42.716 45.9340  51.78545  47.9820  58.0755   76.069    20
##    eigen_PivQR_lapacke 48.567 52.8090  62.94680  62.0250  73.5820   86.894    20
##     eigen_LDLt_lapacke 43.301 47.5430  55.32545  50.4685  64.9515   75.483    20
##      eigen_SVD_lapacke 84.845 87.7715 109.61140 106.9340 127.5605  153.891    20
##  eigen_SymmEig_lapacke 64.366 65.5360  81.86140  71.6800  87.0400  203.629    20
##       eigen_QR_lapacke 46.226 47.5435  56.68590  52.0780  67.2915   79.287    20
##      eigen_LLt_lapacke 47.689 49.8840  56.58340  51.3465  65.2430   74.898    20
##            arma_solve1 91.575 94.7925 114.40970  97.7185 132.9735  177.005    20
##            arma_solve2 28.379 33.2075  85.16750  35.9865  45.6415  970.160    20
##                arma_qr 28.379 35.8400 100.13240  45.2025  49.2985 1192.220    20
##              arma_chol 24.576 28.8190  93.78365  31.1595  36.7175 1239.030    20
##              arma_pinv 44.471 46.2265 113.98540  47.3960  63.9270 1205.970    20
##                    GSL 42.131 43.3010  52.56070  48.2745  59.8310   75.483    20

sessionInfo()
# Reference benchmark runs on windows 10 x64 with MRO-3.4.0,
# RcppGSL 0.3.2, RcppArmadillo 0.7.800.2.0 and RcppEigen 0.3.3.3.0
# and CPU is i7-3770K@4.2GHz.

.Call("RcppBlaze_blaze_version", FALSE, PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_SSE", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_AVX", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_AVX2", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_MIC", PACKAGE = "RcppBlaze")

.Call("RcppBlaze_Blaze_FMA", PACKAGE = "RcppBlaze")
