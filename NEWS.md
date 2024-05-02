# RcppBlaze v1.0.1

* Add function to set/get the number of threads.
* Add function to set/get seed for blaze.
* Improve the performance of `fastLm` by using `CustomVector` and `CustomMatrix`.

# RcppBlaze v1.0.0

* Update Blaze-lib to 3.8.2 since Rcpp supports C++14 / C++17 since 4.2.0.
* Because of the support of C++14, R need to be >= `4.2.0` to install `RcppBlaze`.
* Enable C++11 Threads when compiling `fastLm`.
* Change to use BSD License which is aligned with `blize-lib`.

# RcppBlaze v0.2.2

* Fix building problem on windows.

# RcppBlaze v0.2.1

* Modification of examples.

# RcppBlaze v0.2.0

* Make RcppBlaze link to RcppLAPACKE.
* Add examples to link RcppLAPACKE.

# RcppBlaze v0.1.1

* Fix some minor bugs.

# RcppBlaze v0.1.0

* First release on CRAN.
