#' R and Blaze Integration
#'
#' RcppBlaze construct a bridge between R and Blaze.
#'
#' 'Blaze' is an open-source, high-performance C++ math library for dense and sparse arithmetic.
## Copyright (C) 2010 - 2024 Dirk Eddelbuettel, Romain Francois and Douglas Bates
## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is based on files from RcppArmadillo.
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

#' With its state-of-the-art Smart Expression Template implementation \strong{Blaze} combines the elegance and
#' ease of use of a domain-specific language with HPC-grade performance, making it one of the most
#' intuitive and fastest C++ math libraries available. The \strong{RcppBlaze} package includes the header files
#' from the \strong{Blaze} library with disabling some functionalities related to link to the thread and system
#' libraries which make \strong{RcppBlaze} be a header-only library. Therefore, users do not need to  install
#' \strong{Blaze}.
#'
#' @section Using RcppBlaze:
#' The simplest way to get started is to create a skeleton of a package
#' using \code{RcppBlaze}. This can be done conveniently by the
#' \code{\link{RcppBlaze.package.skeleton}} function.
#'
#' The important steps are
#' \enumerate{
#' \item Include the \samp{RcppBlaze.h} header file, which also includes \samp{blaze/Blaze.h}.
#' \item Import \code{Rcpp}, LinkingTo \code{Rcpp} and \code{RcppBlaze} by adding these lines to the \samp{DESCRIPTION} file:
#' \preformatted{
#'   Imports: Rcpp (>= 1.0.0)
#'   LinkingTo: Rcpp, RcppBlaze
#' }
#' \item Link against the \code{BLAS} and \code{LAPACK} libraries, by adding following two lines in the \samp{Makevars} and \samp{Makevars.win} files:
#' \preformatted{
#'   PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
#' }
#' }
#'
#' Note that if you would like to enable Boost threads support, you need to import \strong{BH} package in your DESCRIPTION.
#' Note that \code{CompressedVector} and \code{CompressedMatrix} only support \code{int}, \code{float} and \code{double} types.
#'
#' @author
#' For RcppBlaze: Ching-Chuan Chen
#' Maintainer: Ching-Chuan Chen <zw12356@gmail.com>
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
#' @name RcppBlaze-package
#' @useDynLib RcppBlaze, .registration = TRUE
#' @importClassesFrom Matrix dgCMatrix lgCMatrix dgRMatrix lgRMatrix
"_PACKAGE"
