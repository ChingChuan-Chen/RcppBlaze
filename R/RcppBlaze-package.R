#' R and Blaze Integration
#'
#' RcppBlaze construct a bridge between R and Blaze.
#'
#' 'Blaze' is an open-source, high-performance C++ math library for dense and sparse arithmetic.
## Copyright (C) 2010 - 2013 Dirk Eddelbuettel, Romain Francois and Douglas Bates
## Copyright (C) 2014        Dirk Eddelbuettel
## Copyright (C) 2017        Chingchuan Chen
##
## This file is based on RcppEigen-package.Rd and RcppArmadillo-package.Rd
## from RcppArmadillo and RcppEigen.
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppBlaze is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppBlaze  If not, see <http://www.gnu.org/licenses/>.

#' With its state-of-the-art Smart Expression Template implementation \strong{Blaze} combines the elegance and
#' ease of use of a domain-specific language with HPC-grade performance, making it one of the most
#' intuitive and fastest C++ math libraries available. The \strong{Blaze} library offers:
#' \itemize{
#'   \item high performance through the integration of BLAS libraries and manually tuned HPC math kernels
#'   \item vectorization by SSE, SSE2, SSE3, SSSE3, SSE4, AVX, AVX2, AVX-512, FMA, and SVML
#'   \item parallel execution by OpenMP, C++11 threads and \strong{Boost} threads
#'      (\strong{Boost} threads is disables in \strong{RcppBlaze})
#'   \item the intuitive and easy to use API of a domain specific language
#'   \item unified arithmetic with dense and sparse vectors and matrices
#'   \item thoroughly tested matrix and vector arithmetic
#'   \item completely portable, high quality C++ source code
#' }
#'
#' The \strong{RcppBlaze} package includes the header files from the \strong{Blaze} library with disabling some
#' functionalities related to link to the thread and system libraries which make \strong{RcppBlaze} be a
#' header-only library. Therefore, users do not need to  install \strong{Blaze} and the dependency \strong{Boost}.
#' \strong{Blaze} is licensed under the New (Revised) BSD license, while \strong{RcppBlaze}
#' (the \strong{Rcpp} bindings/bridge to \strong{Blaze}) is licensed under the GNU GPL version 2 or later,
#' as is the rest of \strong{Rcpp}.
#'
#' Note that since \strong{Blaze} has committed to C++14 which does not used by most \strong{R} users from version 3.0,
#' we will use the version 2.6 of \strong{Blaze} which is C++98 compatible to support the most compilers and system.
#'
#' @section Caution:
#' On the windows x64 with \strong{Rtools33}, \strong{Rtools34}, \code{std::ptrdiff_s} is defined as \code{long long unsigned int},
#' so do not direct wrap the dimention of matrix or the size of vector to output in case the compiling error occurs.
#'
#' @section Using RcppBlaze:
#' The simplest way to get started is to create a skeleton of a package
#' using \code{RcppBlaze}. This can be done conveniently by the
#' \code{\link{RcppBlaze.package.skeleton}}
#' function.
#'
#' The important steps are
#' \enumerate{
#' \item Include the \samp{RcppBlaze.h} header file, which also includes \samp{blaze/Blaze.h}.
#' \item Import \code{Rcpp}, LinkingTo \code{Rcpp}, \code{BH} and \code{RcppBlaze} by adding these lines to the \samp{DESCRIPTION} file:
#' \preformatted{
#'   Imports: Rcpp (>= 0.11.0)
#'   LinkingTo: Rcpp, BH, RcppBlaze
#' }
#' \item Link against the \code{BLAS} and \code{LAPACK} libraries, by adding following two lines in the \samp{Makevars} and
#' \samp{Makevars.win} files:
#' \preformatted{
#'   PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
#' }
#' }
#'
#' @author
#' For RcppBlaze: Chingchuan Chen
#' Maintainer: Chingchuan Chen <zw12356@gmail.com>
#' For blaze: Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff
#'
#' @references
#' \enumerate{
#' \item Blaze project: \url{https://bitbucket.org/blaze-lib/blaze}
#' \item K. Iglberger, G. Hager, J. Treibig, and U. Ruede:
#'       \href{http://epubs.siam.org/sisc/resource/1/sjoce3/v34/i2/pC42_s1}{
#'         Expression Templates Revisited: A Performance Analysis of Current Methodologies.}
#'       SIAM Journal on Scientific Computing, 34(2): C42--C69, 2012
#' \item K. Iglberger, G. Hager, J. Treibig, and U. Ruede:
#'       \href{http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=06266939}{
#'       High Performance Smart Expression Template Math Libraries.}
#'       Proceedings of the 2nd International Workshop on New Algorithms and Programming Models
#'       for the Manycore Era (APMM 2012) at HPCS 2012
#' }
#'
#' @keywords package interface
#' @docType package
#' @name RcppBlaze-package
#' @useDynLib RcppBlaze
#' @importFrom Rcpp evalCpp
#' @importClassesFrom Matrix dgCMatrix lgCMatrix dgRMatrix lgRMatrix
NULL
