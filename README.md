## RcppBlaze

![blaze300x150.jpg](https://bitbucket.org/blaze-lib/blaze/wiki/images/blaze300x150.jpg)

[Blaze](https://bitbucket.org/blaze-lib/blaze) is an open-source, high-performance C++ math library 
for dense and sparse arithmetic. With its state-of-the-art *Smart Expression Template* implementation 
**Blaze** combines the elegance and   ease of use of a domain-specific language with HPC-grade performance, 
making it one of the most  intuitive and fastest C++ math libraries available. The **RcppBlaze** package includes 
the header files from the **Blaze** library with disabling some functionalities related to link to the thread and 
system libraries which make **RcppBlaze** be a header-only library. Therefore, users do not need to  install 
**Blaze**. **Blaze** is licensed under the New (Revised) BSD license, while **RcppBlaze** 
(the 'Rcpp' bindings/bridge to **Blaze**) is licensed under the GNU GPL version 2 or later, as is the rest of **Rcpp**. 

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

### Authors

Ching-Chuan Chen, Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff

### License

GPL (>= 2)
