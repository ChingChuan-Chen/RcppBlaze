## RcppBlaze

[![Build Status](https://travis-ci.org/ChingChuan-Chen/RcppBlaze.svg)](https://travis-ci.org/ChingChuan-Chen/RcppBlaze) [![Build status](https://ci.appveyor.com/api/projects/status/ip1s0avhreksvgd5?svg=true)](https://ci.appveyor.com/project/ChingChuan-Chen/rcppblaze) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

### Overview

![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

[Blaze](https://bitbucket.org/blaze-lib/blaze) is an open-source, high-performance C++ math library 
for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation 
**Blaze** combines the elegance and   ease of use of a domain-specific language with HPC-grade performance, 
making it one of the most  intuitive and fastest C++ math libraries available. The **Blaze** library offers:

   * high performance through the integration of BLAS libraries and manually tuned HPC math kernels
   * vectorization by SSE, SSE2, SSE3, SSSE3, SSE4, AVX, AVX2, AVX-512, FMA, and SVML
   * parallel execution by OpenMP, C++11 threads and **Boost** threads (**Boost** threads is disables in **RcppBlaze**)
   * the intuitive and easy to use API of a domain specific language
   * unified arithmetic with dense and sparse vectors and matrices
   * thoroughly tested matrix and vector arithmetic
   * completely portable, high quality C++ source code
   
The **RcppBlaze** package includes the header files from the **Blaze** library with disabling some
functionalities related to link to the thread and system libraries which make **RcppBlaze** be a 
header-only library. Therefore, users do not need to  install **Blaze** and the dependency **Boost**. 
**Blaze** is licensed under the New (Revised) BSD license, while **RcppBlaze**
(the 'Rcpp' bindings/bridge to **Blaze**) is licensed under the GNU GPL version 2 or later, 
as is the rest of **Rcpp**. 

Note that since \strong{Blaze} has committed to C++14 which does not used by most R users from version 3.0, 
we will use the version 2.6 of **Blaze** which is C++98 compatible to support the most compilers and system.

### Installation

You can install:

* the stable version from CRAN with

    ```R
    install.packages("RcppBlaze")
    ```

* the latest development version from github with

    ```R
    install.packages("devtools")
    devtools::install_github("ChingChuan-Chen/RcppBlaze")
    ```

If you encounter a bug, please file a reproducible example on [github](https://github.com/ChingChuan-Chen/RcppBlaze/issues).

### Authors

Chingchuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

GPL (>= 2)
