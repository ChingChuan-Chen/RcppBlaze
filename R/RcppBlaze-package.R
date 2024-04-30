## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the 3-Clause BSD License. You should have received
## a copy of 3-Clause BSD License along with RcppBlaze.
## If not, see https://opensource.org/license/BSD-3-Clause.

#' RcppBlaze - 'Rcpp' Integration for the 'Blaze' High-Performance 'C++' Math Library
#'
#' \strong{RcppBlaze} constructs a bridge between \strong{R} and \strong{Blaze}.
#'
#' \strong{Blaze} is an open-source, high-performance \strong{C++} math library for dense and sparse arithmetic.
#' With its state-of-the-art Smart Expression Template implementation \strong{Blaze} combines the elegance and
#' ease of use of a domain-specific language with HPC-grade performance, making it one of the most
#' intuitive and fastest \strong{C++} math libraries available. The \strong{RcppBlaze} package includes the header files
#' from the \strong{Blaze} library with disabling some functionalities related to link to the thread and system
#' libraries which make \strong{RcppBlaze} be a header-only library. Therefore, users do not need to  install
#' \strong{Blaze}.
#'
#' @section Using \strong{RcppBlaze}:
#' To use \strong{RcppBlaze} in your package, there are some important steps:
#' \enumerate{
#' \item Include the \samp{RcppBlaze.h} header file, which also includes \samp{blaze/Blaze.h}.
#' \item Import \code{Rcpp}, LinkingTo \code{Rcpp} and \code{RcppBlaze} by adding these lines to the \samp{DESCRIPTION} file:
#' \preformatted{
#'   Imports: Rcpp (>= 1.0.0)
#'   LinkingTo: Rcpp, RcppBlaze
#' }
#' \item Link against the \code{BLAS} and \code{LAPACK} libraries, by adding following two lines in the \samp{Makevars} and \samp{Makevars.win} files:
#' \preformatted{
#'   PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)
#'   PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CXXFLAGS)
#' }
#' \item Since there are conflicted definitions between \strong{R} and \strong{blaze} which is \code{TRUE} and \code{FALSE}.
#'  You have to write the initializing function for \strong{C/C++} code which the function is named after \code{R_init_YourPackageName}
#'  You can refer to our another package, \url{https://github.com/ChingChuan-Chen/RcppLbfgsBlaze} for example.
#' }
#'
#' @section Notes:
#' \enumerate{
#' \item If you would like to enable Boost threads support, you need to import \strong{BH} package in your DESCRIPTION.
#' \item \code{CompressedVector} and \code{CompressedMatrix} only support \code{int}, \code{float} and \code{double} types.
#' }
#'
#' @author
#' For RcppBlaze: Ching-Chuan Chen
#' Maintainer: Ching-Chuan Chen <zw12356@gmail.com>
#' For blaze: Klaus Iglberger, Georg Hager, Christian Godenschwager, Tobias Scharpff
#'
#' @references
#' \enumerate{
#' \item Blaze project: \url{https://bitbucket.org/blaze-lib/blaze}.
#' \item K. Iglberger, G. Hager, J. Treibig, and U. Ruede:
#'       Expression Templates Revisited: A Performance Analysis of Current Methodologies.
#'       SIAM Journal on Scientific Computing, 34(2): C42--C69, 2012, \doi{10.1137/110830125}.
#' \item K. Iglberger, G. Hager, J. Treibig, and U. Ruede,
#'       High Performance Smart Expression Template Math Libraries.
#'       Proceedings of the 2nd International Workshop on New Algorithms and Programming Models
#'       for the Manycore Era (APMM 2012) at HPCS 2012, \doi{10.1109/HPCSim.2012.6266939}.
#' }
#'
#' @keywords package interface
#' @name RcppBlaze-package
#' @useDynLib RcppBlaze, .registration = TRUE
#' @importClassesFrom Matrix dgCMatrix lgCMatrix dgRMatrix lgRMatrix
"_PACKAGE"
