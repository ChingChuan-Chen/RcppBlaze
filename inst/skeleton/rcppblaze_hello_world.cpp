// Copyright (C) 2010 - 2013 Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C) 2014        Dirk Eddelbuettel
// Copyright (C)  2017 - 2024  Ching-Chuan Chen
//
// This file is based rcppeigen_hello_world.cpp and
// rcpparma_hello_world.cpp from RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.

// we only include RcppBlaze.h which pulls Rcpp.h in for us
#include "RcppBlaze.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppBlaze so that the build process will know what to do
//
// [[Rcpp::depends(RcppBlaze)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
blaze::DynamicMatrix<int> rcppblaze_hello_world() {
  blaze::DiagonalMatrix< blaze::DynamicMatrix<int> > m1( 3 ), m2( 3 );
  for( size_t i=0UL; i<3UL; ++i) {
    m1(i,i) = 1;
    m2(i,i) = 1;
  }
  return m1 + 3 * (m1 + m2);
}

// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
blaze::DynamicMatrix<double> rcppblaze_outerproduct(const blaze::DynamicVector<double>& x) {
  blaze::DynamicMatrix<double> m = x * blaze::trans(x);  // Outer product between two vectors
  return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcppblaze_innerproduct(const blaze::DynamicVector<double>& x) {
  return (x, x);
}

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcppblaze_bothproducts(const blaze::DynamicVector<double>& x) {
  blaze::DynamicMatrix<double> op = x * blaze::trans(x);
  return Rcpp::List::create(
    Rcpp::Named("outer") = op,
    Rcpp::Named("inner") = (x, x)
  );
}
