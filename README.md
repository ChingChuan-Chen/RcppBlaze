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

    ```R
    install.packages("RcppBlaze")
    ```

* the latest development version from github with

    ```R
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

    ```R
    source(system.file("examples", "lmBenchmark.R", package = "RcppBlaze"))
    # lm benchmark for n = 10000 and p = 100: nrep = 20
    # Unit: milliseconds
    #               expr      min        lq       mean    median        uq      max neval
    #             lm.fit  28.4741  30.13910  34.926360  30.84910  32.18480 108.9459    20
    #           blaze_qr  64.5394  68.55650  69.941730  69.87310  72.03110  74.3638    20
    #         blaze_ldlt   2.6675   3.02795   3.795505   3.61120   4.17630   6.2245    20
    #          blaze_llt   2.7603   3.22650   3.910560   3.70695   4.32075   6.2395    20
    #        eigen_PivQR  17.1567  17.72335  18.982210  18.94920  19.58690  24.6407    20
    #         eigen_LDLt   2.1242   2.26725   2.451755   2.41675   2.52845   3.0334    20
    #          eigen_SVD 122.6997 126.32880 129.656900 129.21315 133.09995 139.8186    20
    #      eigen_SymmEig   7.0235   7.51910   7.896200   7.86275   8.24215   8.9439    20
    #       lapack_GESDD  96.4349 101.48350 103.219585 102.96085 105.07070 110.5083    20
    #           eigen_QR  12.2910  12.88645  14.254070  13.89960  14.67575  24.2232    20
    #          eigen_LLt   2.1843   2.26675   2.435735   2.38790   2.64215   2.8397    20
    #        arma_fastLm  53.2420  55.59120  56.801435  56.41610  57.72490  64.7978    20
    #  arma_direct_solve  18.4159  19.63645  20.332250  20.04015  20.80255  23.8401    20
    #            arma_qr  62.8046  64.96840  67.617250  66.71465  69.89655  76.6999    20
    #          arma_chol  18.7052  19.50140  20.039115  19.97345  20.42390  22.2089    20
    #          arma_pinv  19.5845  20.20290  21.070095  21.01790  21.57585  24.5399    20
    #                GSL 225.8466 278.33925 294.642770 300.72225 315.53745 337.1306    20
    ```

Above results are run on my desktop (i9-13900K, DDR5-4000 128GB).

### Authors

Ching-Chuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

BSD-3 License
