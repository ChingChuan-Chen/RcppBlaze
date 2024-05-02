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
    #             lm.fit  28.1905  29.09360  33.631160  30.26040  31.51810  98.2644    20
    #           blaze_qr  65.1809  66.35890  69.649720  67.53340  68.74100 109.3559    20
    #         blaze_ldlt   2.9600   3.24410   4.057645   3.86770   4.62575   5.5541    20
    #          blaze_llt   2.9146   3.50550   4.202165   4.14505   5.03215   5.6853    20
    #        eigen_PivQR  17.1320  17.52165  18.378385  17.94195  18.90675  22.4894    20
    #         eigen_LDLt   2.1492   2.23390   2.399370   2.35980   2.47280   3.4185    20
    #          eigen_SVD 120.5626 122.29585 129.181575 124.52630 126.36725 199.2691    20
    #      eigen_SymmEig   6.9530   7.15720   7.666565   7.62780   7.89635   9.1642    20
    #       lapack_GESDD  95.3666  97.76600  99.878040  99.61240 100.62335 113.1822    20
    #           eigen_QR  12.1408  12.61675  13.241620  13.27330  13.86170  14.4592    20
    #          eigen_LLt   2.0713   2.16255   2.456665   2.29280   2.41630   5.2753    20
    #        arma_fastLm  52.3471  53.28430  56.630075  54.01670  55.92150  95.4131    20
    #  arma_direct_solve  18.6326  18.93075  19.547790  19.14285  19.73080  24.3872    20
    #            arma_qr  63.6403  65.18815  66.846170  66.70980  67.29430  77.9830    20
    #          arma_chol  18.5586  19.32225  20.054260  19.73960  20.61855  22.4312    20
    #          arma_pinv  19.6419  20.04995  20.525605  20.51070  20.78480  22.0781    20
    #                GSL 231.2039 268.38020 291.832665 299.65985 316.51120 352.4575    20
    ```

Above results are run on my desktop (i9-13900K, DDR5-4000 128GB).

### Authors

Ching-Chuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

BSD-3 License
