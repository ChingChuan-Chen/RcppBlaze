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

exprs$blaze_qr <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', mm, y, 0L))
exprs$blaze_ldlt <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', mm, y, 1L))
exprs$blaze_llt <- expression(.Call('_RcppBlaze_fastLmPure', PACKAGE = 'RcppBlaze', mm, y, 2L))

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

  if (suppressMessages(require("RcppLAPACKE", character = TRUE, quietly = TRUE)) &&
      packageVersion("RcppEigen") >= "0.3.3.3.0") {
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
##                   expr         min          lq       mean     median         uq        max neval
##                 lm.fit  161.376588  182.397185  264.96363  228.84149  335.04137  499.73649    20
##               blaze_qr  127.149891  143.475175  284.87366  188.23037  436.66301  750.07794    20
##             blaze_ldlt   41.078946   44.108944   48.61167   46.37957   51.61782   66.14742    20
##              blaze_llt   40.642524   42.656292   46.50010   45.68276   49.14622   61.81600    20
##         blaze_qr_cblas  132.741600  149.980077  228.39954  180.38604  268.08239  620.94724    20
##       blaze_ldlt_cblas   41.057845   42.000726   45.44590   43.55681   47.75977   60.06425    20
##        arma_chol_cblas   40.884076   41.324976   43.76919   43.21980   45.95689   48.13963    20
##            eigen_PivQR  135.794738  137.935854  147.88412  140.05725  147.18336  255.20435    20
##             eigen_LDLt   28.666661   29.298621   30.46676   29.76745   30.64507   36.32741    20
##              eigen_SVD 1071.651124 1207.986435 1314.85447 1304.42064 1414.64429 1641.33690    20
##          eigen_SymmEig   71.522409   73.605271   75.92102   74.84183   75.66863   90.01379    20
##           lapack_GESDD  188.888941  230.760178  294.28531  246.34024  377.51085  520.11067    20
##               eigen_QR  128.969048  134.426637  141.10472  138.30824  145.77743  166.68486    20
##              eigen_LLt   28.127403   29.007185   29.94602   29.95371   30.37272   33.52035    20
##    eigen_PivQR_lapacke  142.557837  158.791849  186.14268  180.58448  195.75288  291.02143    20
##     eigen_LDLt_lapacke    8.586206    9.361428   14.50333   15.10735   18.51547   21.61671    20
##      eigen_SVD_lapacke  149.135728  274.490508  441.21498  343.54774  602.71209  929.08686    20
##  eigen_SymmEig_lapacke   19.607602   23.866080   28.16363   26.05714   29.90145   47.68422    20
##       eigen_QR_lapacke   42.896236   56.873169   98.34746   72.52742  101.51413  348.19020    20
##      eigen_LLt_lapacke    8.742399   10.053396   14.15510   12.62205   16.87693   24.99200    20
##            arma_solve1   55.563868   73.501424   99.98713   83.75828  120.96035  202.57762    20
##            arma_solve2    8.489917   11.595660   16.08504   15.30364   21.58805   24.84006    20
##                arma_qr  101.442060  114.226044  152.15244  126.20073  149.53254  377.03877    20
##              arma_chol    8.313652    9.254720   13.75527   12.46967   15.30539   33.37284    20
##              arma_pinv    8.883698   10.186594   16.35474   13.25469   18.38365   38.40573    20
##                    GSL 3734.103710 3822.975717 3911.08668 3881.08361 3981.66668 4132.22756    20

print(do_bench(n = 1.6e4L, p = 250L))
# Reference benchmark:
## Unit: milliseconds
##                   expr       min        lq       mean     median         uq        max neval
##                 lm.fit 829.90126 957.59911 1397.10646 1057.85109 1253.27194 5363.31948    20
##               blaze_qr 253.01198 270.92673  321.24955  295.65067  361.10560  588.83688    20
##             blaze_ldlt  69.22110  85.44591   96.71150   96.31355  105.00747  127.62336    20
##              blaze_llt  74.44316  79.90160   96.09635   90.90350  108.87658  146.38436    20
##         blaze_qr_cblas 258.59032 286.38827  356.62706  305.19216  400.68590  600.82498    20
##       blaze_ldlt_cblas 190.96416 193.32769  203.48333  200.03642  209.84329  234.47938    20
##        arma_chol_cblas 189.40559 192.76065  207.74746  198.38357  219.36541  260.79429    20
##            eigen_PivQR 631.70973 645.84079  669.36712  652.79332  693.77426  751.21490    20
##             eigen_LDLt 118.18805 120.04018  122.63647  121.12950  123.35789  143.33302    20
##          eigen_SymmEig 365.65014 368.96383  374.99892  371.77591  378.55913  393.82936    20
##           lapack_GESDD 508.51234 550.99333  687.21797  693.85308  784.60552  948.71844    20
##               eigen_QR 298.38563 300.79958  309.38828  304.21814  315.70761  340.39650    20
##              eigen_LLt 119.05636 120.78376  123.18109  122.93745  125.08235  132.07293    20
##    eigen_PivQR_lapacke 393.01754 519.32013  653.93112  664.55172  759.79825  971.30089    20
##     eigen_LDLt_lapacke  18.94047  22.12653   27.78196   28.47185   33.25617   41.34771    20
##      eigen_SVD_lapacke 351.36247 428.99136  535.31808  500.43238  642.28390  830.99344    20
##  eigen_SymmEig_lapacke  57.41224  66.90229   98.70740   78.86624  115.23242  239.09326    20
##       eigen_QR_lapacke  78.41580 102.26411  171.02293  138.78060  224.60316  424.00678    20
##      eigen_LLt_lapacke  20.10112  21.79846   32.75467   26.99539   34.94872   96.45657    20
##            arma_solve1 123.29151 144.38238  239.41814  168.12967  304.94795  604.82268    20
##            arma_solve2  22.12124  26.90400   39.28202   32.53313   46.56711   85.73000    20
##                arma_qr 234.41134 283.07752  342.38603  312.47686  372.67618  596.41297    20
##              arma_chol  19.46355  20.06615   35.80455   25.13173   41.89529   89.81933    20
##              arma_pinv  32.43387  36.66914   61.95252   42.58220   71.38084  157.53695    20

print(do_bench(n = 50, p = 10))
# Reference benchmark:
## Unit: microseconds
##                   expr     min       lq      mean   median       uq      max neval
##                 lm.fit  85.999 100.4835 135.16285 122.1285 162.5175  249.065    20
##               blaze_qr  25.949  28.7185  34.89855  33.1415  39.7905   48.551    20
##             blaze_ldlt  11.381  17.8275  23.13805  20.1745  23.1080   55.072    20
##              blaze_llt  17.099  19.8425  22.49565  20.6035  24.0785   34.287    20
##         blaze_qr_cblas  21.231  26.7300  33.63270  28.4010  34.1000  102.464    20
##       blaze_ldlt_cblas  16.630  18.8155  25.16870  23.2545  24.5345   52.235    20
##        arma_chol_cblas  14.656  19.5390  23.58990  21.6075  25.8190   44.851    20
##            eigen_PivQR  51.377  53.6805  67.74720  65.6310  76.3700   99.920    20
##             eigen_LDLt  42.607  47.6520  72.33745  62.3100  85.0135  223.165    20
##              eigen_SVD  93.759  97.4790 120.69455 112.3050 141.0705  165.538    20
##          eigen_SymmEig  53.618  61.8410  75.90005  74.4520  82.4430  124.203    20
##           lapack_GESDD 105.814 114.4095 143.18730 137.9170 163.0690  238.920    20
##               eigen_QR  44.420  49.4495  69.37080  61.5060  76.9500  143.783    20
##              eigen_LLt  44.619  47.4780  58.46490  54.3350  70.2970   85.673    20
##    eigen_PivQR_lapacke  62.366  67.0140  79.48070  79.3535  86.6435  111.646    20
##     eigen_LDLt_lapacke  51.390  53.6665  69.72835  64.7315  72.6190  145.433    20
##      eigen_SVD_lapacke 105.712 107.2895 128.61790 119.7340 144.8710  180.529    20
##  eigen_SymmEig_lapacke  73.608  77.4120  99.08870  87.3300 103.4105  218.399    20
##       eigen_QR_lapacke  52.919  55.2880  70.23635  66.0800  76.5300  128.366    20
##      eigen_LLt_lapacke  52.336  58.6425  75.01505  72.4995  86.3405  127.881    20
##            arma_solve1 173.182 184.3575 898.93265 233.0770 865.4830 4756.314    20
##            arma_solve2  38.007  43.6945  57.19055  54.6490  61.7275  113.098    20
##                arma_qr  31.683  40.3385  53.37020  45.3865  57.5395  145.428    20
##              arma_chol  35.992  38.9905  50.66730  48.9505  54.5160   90.422    20
##              arma_pinv  65.466  67.5415  82.59380  81.2200  93.2570  119.064    20
##                    GSL  63.371  70.8065  91.29635  88.8325  95.1760  193.896    20

sessionInfo()
# Reference benchmarks are run on CentOS 7 x86_64 with MRO-3.4.0,
# RcppGSL 0.3.2, RcppArmadillo 0.7.900.2.0 and RcppEigen 0.3.3.3.0.
# The machine is equiped with CPU is i7-3770K@4.2GHz.

.Call("_RcppBlaze_blaze_version", FALSE, PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_SSE", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_AVX", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_AVX2", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_MIC", PACKAGE = "RcppBlaze")

.Call("_RcppBlaze_Blaze_FMA", PACKAGE = "RcppBlaze")
