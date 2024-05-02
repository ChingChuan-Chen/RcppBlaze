## RcppBlaze

![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

[Blaze](https://bitbucket.org/blaze-lib/blaze) is an open-source, high-performance **C++** math library 
for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation 
**Blaze** combines the elegance and ease of use of a domain-specific language with HPC-grade performance, 
making it one of the most intuitive and fastest **C++** math libraries available. The **RcppBlaze** package includes 
the header files from the **Blaze** library with disabling some functionalities related to link to the thread and 
system libraries which make **RcppBlaze** be a header-only library. Therefore, users do not need to  install 
**Blaze**.

### Installation

You can install:

* the stable version from CRAN with

  ``` r
  install.packages("RcppBlaze")
  ```

* the latest development version from github with

  ``` r
  install.packages("remotes")
  remotes::install_github("ChingChuan-Chen/RcppBlaze")
  ```

If you encounter a bug, please file a reproducible example on [github](https://github.com/ChingChuan-Chen/RcppBlaze/issues).

### Note

#### Sparse Matrix only supports three data types
`CompressedVector` and `CompressedMatrix` only support `int`, `float` and `double` types.
Logical part is not supported because there is no enough resource for it (helps are welcome).
Since `Matrix` only provides `l*[CTR]Matrix` and `d*[CTR]Matrix` and does not support `z*[CTR]Matrix`.

#### CustomVector and CustomMatrix cannot directly be converted from R object

Since `CustomVector` and `CustomMatrix` are mapping memory which user should manage by himself.
Therefore, we provide another function to help user do convertion from R objects.
We provide `RcppBlaze::copyToCustomVector` and `RcppBlaze::copyToCustomMatrix` to make easier on data copy. 
Below is the example code:

``` c++
// For CustomVector<int>
using iCustomVectorUU = blaze::CustomVector<int, blaze::unaligned, blaze::unpadded>;
using iCustomVectorAP = blaze::CustomVector<int, blaze::aligned, blaze::padded>;

// initialize R IntegerVector
Rcpp::IntegerVector intVec = Rcpp::IntegerVector::create(-2, -1, 0, 1, 2 );
size_t int_vec_size = Rf_xlength(intVec);
size_t intVecPaddedSize = blaze::nextMultiple<size_t>(int_vec_size, blaze::SIMDTrait<int>::size);

// unaligned & unpadded CustomVector
std::unique_ptr<int[], blaze::ArrayDelete> data_unpadded(new int[int_vec_size]);
iCustomVectorUU cv_uu_int(data_unpadded.get(), int_vec_size);
RcppBlaze::copyToCustomVector(intVec, cv_uu_int);

// aligned & padded CustomVector
std::unique_ptr<int[], blaze::Deallocate> data_padded(blaze::allocate<int>(intVecPaddedSize));
iCustomVectorAP cv_ap_int(data_padded.get(), int_vec_size, intVecPaddedSize);
RcppBlaze::copyToCustomVector(intVec, cv_ap_int);

// For CustomMatrix<int>
using iCustomMatrixUU = blaze::CustomMatrix<int, blaze::unaligned, blaze::unpadded, blaze::columnMajor>;
using iCustomMatrixAP = blaze::CustomMatrix<int, blaze::aligned, blaze::padded, blaze::columnMajor>;
using iCustomMatrixAP_RM = blaze::CustomMatrix<int, blaze::aligned, blaze::padded, blaze::rowMajor>;

// initialize R IntegerMatrix
Rcpp::IntegerMatrix intMat(2, 3);
std::fill(intMat.begin(), intMat.end(), 8);
intMat(0, 1) = 5;
intMat(1, 2) = 4;
Rcpp::Shield<SEXP> intMatDimsSexp(Rf_getAttrib(intMat, R_DimSymbol));
int* intMatDims = INTEGER(intMatDimsSexp);
size_t m = (size_t) intMatDims[0], n = (size_t) intMatDims[1];

// column-major parameters
size_t intSimdSize = blaze::SIMDTrait<int>::size;
size_t intMatPaddedRows = blaze::nextMultiple<size_t>(m, intSimdSize);

// unaligned & unpadded column-major CustomMatrix
std::unique_ptr<int[], blaze::ArrayDelete> data_unpadded(new int[m*n]);
iCustomMatrixUU cm_uu_int(data_unpadded.get(), m, n);
RcppBlaze::copyToCustomMatrix(intMat, cm_uu_int);

// aligned & padded column-major CustomMatrix
std::unique_ptr<int[], blaze::Deallocate> data_padded(blaze::allocate<int>(intMatPaddedRows * n));
iCustomMatrixAP cm_ap_int(data_padded.get(), m, n, intMatPaddedRows);
RcppBlaze::copyToCustomMatrix(intMat, cm_ap_int);

// row-major parameters
size_t intMatPaddedCols = blaze::nextMultiple<size_t>(n, intSimdSize);

// aligned & padded row-major CustomMatrix
std::unique_ptr<int[], blaze::Deallocate> data_rm_padded(blaze::allocate<int>(m * intMatPaddedCols));
iCustomMatrixAP_RM cm_ap_rm_int(data_rm_padded.get(), m, n, intMatPaddedCols);
RcppBlaze::copyToCustomMatrix(intMat, cm_ap_rm_int);
```

### Linear Model Fitting Benchmark 

You can refer to the file [lmBenchmark.R](./inst/examples/lmBenchmark.R) to find the code.
Below code and corresponding results show that `RcppBlaze` have better performance than `RcppArmadillo` and `RcppGSL`.
However, `RcppEigen` can provide more efficient algorithms (`LDLt` and `LLt`) for linear model fitting (about 2.4 times faster).

``` r
source(system.file("examples", "lmBenchmark.R", package = "RcppBlaze"))
# lm benchmark for n = 10000 and p = 100: nrep = 20
# Unit: milliseconds
#               expr      min        lq       mean    median        uq      max neval
#             lm.fit  28.9307  30.76115  34.718960  31.32925  32.71110  94.3755    20
#           blaze_qr  66.7372  67.63715  68.735135  68.52440  69.61730  71.2159    20
#         blaze_ldlt   2.7382   3.21195   3.878950   3.48345   3.94490   9.0361    20
#          blaze_llt   2.8695   3.26960   3.807020   3.61675   3.91585   6.9932    20
#           blaze_lu   2.8755   3.12790   3.573530   3.45100   3.78790   5.0707    20
#        eigen_PivQR  16.8238  17.60885  18.730360  19.01100  19.52710  20.0469    20
#         eigen_LDLt   1.9877   2.15395   2.342755   2.28880   2.45690   3.1641    20
#          eigen_SVD 122.9787 127.08120 130.421575 129.57880 133.20360 141.0751    20
#      eigen_SymmEig   7.1041   7.59225   7.911580   7.94090   8.21375   8.9926    20
#       lapack_GESDD  97.5447  99.89155 103.044495 102.50625 104.55925 117.6772    20
#           eigen_QR  12.5883  13.26180  13.731200  13.67780  14.13895  15.0327    20
#          eigen_LLt   2.0867   2.24700   2.384330   2.40915   2.50055   2.6669    20
#        arma_fastLm  52.9779  55.86840  58.040880  57.72055  58.99520  71.6400    20
#  arma_direct_solve  18.6251  19.64045  19.992785  19.96915  20.28160  21.4155    20
#            arma_qr  63.1927  66.41025  68.602565  67.99540  69.83185  76.9579    20
#          arma_chol  18.6167  19.47850  20.352435  19.99895  20.83785  24.0279    20
#          arma_pinv  19.4529  20.48360  20.834835  20.89880  21.03085  22.9177    20
#                GSL 281.8348 294.52890 310.380445 308.12135 325.44070 345.3462    20
```

Above results are run on my desktop (i9-13900K, DDR5-4000 128GB).

### Authors

Ching-Chuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

BSD-3 License
