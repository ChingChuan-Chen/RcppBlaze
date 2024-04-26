## RcppBlaze

![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

[Blaze](https://bitbucket.org/blaze-lib/blaze) is an open-source, high-performance C++ math library 
for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation 
**Blaze** combines the elegance and   ease of use of a domain-specific language with HPC-grade performance, 
making it one of the most  intuitive and fastest C++ math libraries available. The **RcppBlaze** package includes 
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

`CompressedVector` and `CompressedMatrix` only support `int`, `float` and `double` types.
Logical part is not supported because there is no enough resource for it (helps are welcome).
Since `Matrix` only provides `l*[CTR]Matrix` and `d*[CTR]Matrix` and does not support `z*[CTR]Matrix`.

### Known Issues

Converting R object to `CustomVector` or `CustomMatrix` (`Rcpp::as`) is unstable. You should take your own risk to do so.
Returning `CustomVector` or `CustomMatrix` from `C++` (`Rcpp::wrap`) is totally okay.
You can `#define RCPPBLAZE_ENABLE_WARN false` to disable the related warnings.

### Linear Model Fitting Benchmark 

You can refer to the file [lmBenchmark.R](./inst/examples/lmBenchmark.R) to find the code.
Below code and corresponding results show that `RcppBlaze` have better performance than `RcppArmadillo` and `RcppGSL`.
However, `RcppEigen` can provide more efficient algorithms (`LDLt` and `LLt`) for linear model fitting (about 2.4 times faster).

    ```R
    source(system.file("examples", "lmBenchmark.R", package = "RcppBlaze"))
    # lm benchmark for n = 10000 and p = 100: nrep = 20
    # Unit: milliseconds
    #               expr      min        lq       mean    median        uq       max neval
    #             lm.fit  28.2638  29.26775  30.436995  29.74145  30.03235   44.3607    20
    #           blaze_qr  62.7292  67.04560  69.012075  68.70860  70.74640   76.1922    20
    #         blaze_ldlt   4.6197   5.17470   5.956825   5.40375   6.06695    8.8686    20
    #          blaze_llt   4.4422   4.89535   5.804225   5.33140   6.30010    9.1514    20
    #        eigen_PivQR  16.4457  16.70190  17.578635  17.51210  18.63425   18.8380    20
    #         eigen_LDLt   2.0512   2.20940   2.434285   2.51605   2.58100    2.8422    20
    #          eigen_SVD 120.8145 124.06050 126.923975 127.17165 129.10250  133.2594    20
    #      eigen_SymmEig   6.9380   7.25805   7.703965   7.61005   8.08250    8.7749    20
    #       lapack_GESDD  95.6177 100.85300 102.778040 101.75035 103.43840  125.5809    20
    #           eigen_QR  11.9052  13.45050  13.875040  13.96875  14.23940   16.1556    20
    #          eigen_LLt   2.1273   2.30230   2.464870   2.46400   2.62795    2.9289    20
    #        arma_fastLm  53.8013  54.67790  55.902710  55.28985  56.73145   62.8755    20
    #  arma_direct_solve  18.5666  19.44685  20.112800  19.75030  20.35425   23.9168    20
    #            arma_qr  63.5199  64.82850  66.334945  66.19005  67.58970   69.0730    20
    #          arma_chol  18.5555  19.21530  19.711100  19.64705  20.29210   21.0759    20
    #          arma_pinv  19.1042  20.75225  21.186465  20.97250  21.34205   26.8389    20
    #                GSL 221.8449 274.08050 372.676835 287.69530 302.15075 2117.0106    20
    ```

Above results are run on my desktop (i9-13900K, DDR5-4000 128GB).

### Authors

Ching-Chuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

GPL (>= 2)
