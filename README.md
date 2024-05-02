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
    #             lm.fit  27.7852  29.69620  31.400790  31.03340  32.85220  35.4492    20
    #           blaze_qr  64.9408  66.97885  71.221760  70.05315  72.66965  92.3066    20
    #         blaze_ldlt   2.8987   3.14965   3.919910   3.42110   4.20820   8.0185    20
    #          blaze_llt   2.8607   3.16585   3.791535   3.51840   4.19255   6.5455    20
    #        eigen_PivQR  16.8208  17.91055  18.795780  18.56305  19.79840  20.9465    20
    #         eigen_LDLt   2.1944   2.35235   2.598930   2.47610   2.67960   4.6175    20
    #          eigen_SVD 123.3384 128.45530 131.244360 130.01270 133.88340 145.0337    20
    #      eigen_SymmEig   7.1991   7.57645   8.034975   8.01120   8.40675   9.0375    20
    #       lapack_GESDD  95.6156 100.07855 103.519825 102.68900 105.67655 115.6284    20
    #           eigen_QR  12.4275  12.99810  13.980900  13.74980  14.78075  17.1052    20
    #          eigen_LLt   2.1593   2.30635   2.497980   2.50155   2.64545   3.0283    20
    #        arma_fastLm  52.4504  55.20310  56.404615  56.62550  57.76660  60.1031    20
    #  arma_direct_solve  19.1991  19.55165  20.260825  20.04505  20.66995  22.1740    20
    #            arma_qr  61.7838  66.77460  68.316745  68.01125  69.96160  77.8573    20
    #          arma_chol  19.0139  19.54485  20.462105  20.27165  21.00860  23.3993    20
    #          arma_pinv  20.2550  20.90235  21.836590  21.30565  22.76760  26.1785    20
    #                GSL 250.6684 310.45435 324.791035 322.07315 337.41335 409.7735    20
    ```

Above results are run on my desktop (i9-13900K, DDR5-4000 128GB).

### Authors

Ching-Chuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

BSD-3 License
