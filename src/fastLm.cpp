// Copyright (C)  2010 - 2024  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2017 - 2024  Ching-Chuan Chen
//
// This file is based on files from RcppArmadillo.
//
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.

#include <RcppBlaze.h>
using Rcpp::_;

enum {QRSolverType = 0, LDLTSolverType, LLTSolverType};

// [[Rcpp::export]]
Rcpp::List testAs1(Rcpp::NumericMatrix x) {
/*
  typedef typename blaze::CustomMatrix<double, blaze::unaligned, blaze::unpadded, blaze::columnMajor> dCustomMatrixUU;
  dCustomMatrixUU w = Rcpp::as<dCustomMatrixUU>(x);
  Rcpp::Rcout << w << std::endl;

  typedef typename blaze::CustomMatrix<double, blaze::unaligned, blaze::unpadded, blaze::rowMajor> dCustomMatrixUU2;
  dCustomMatrixUU2 y = Rcpp::as<dCustomMatrixUU2>(x);
  Rcpp::Rcout << y << std::endl;
*/

  Rcpp::Rcout << "Column-Major:" << std::endl;
  size_t m = x.rows(), n = x.cols();
  size_t paddedSize = blaze::nextMultiple<size_t>(n, blaze::SIMDTrait<double>::size);
  size_t matSize = (SO == blaze::rowMajor)?(paddedSize*n):(m*paddedSize);
  Rcpp::Rcout << "m: " << m << ", n: " << n << ", paddedSize: " << paddedSize << "matSize: " << matSize << std::endl;
  Rcpp::Rcout << "Row-Major:" << std::endl;
  size_t m2 = x.rows(), n2 = x.cols();
  size_t paddedSize2 = blaze::nextMultiple<size_t>(m, blaze::SIMDTrait<double>::size);
  size_t matSize2 = (SO == blaze::rowMajor)?(paddedSize*n):(m*paddedSize);
  Rcpp::Rcout << "m: " << m2 << ", n: " << n2 << ", paddedSize: " << paddedSize2 << "matSize: " << matSize2 << std::endl;

  typedef typename blaze::CustomMatrix<double, blaze::aligned, blaze::unpadded, blaze::columnMajor> dCustomMatrixAP;

  dCustomMatrixAP z = Rcpp::as<dCustomMatrixAP>(x);
  Rcpp::Rcout << z << std::endl;

  typedef typename blaze::CustomMatrix<double, blaze::aligned, blaze::unpadded, blaze::rowMajor> dCustomMatrixAP2;

  dCustomMatrixAP2 t = Rcpp::as<dCustomMatrixAP2>(x);
  Rcpp::Rcout << t << std::endl;

  return Rcpp::List::create(
    _["test"] = true
  );
}

// [[Rcpp::export]]
Rcpp::List testWrap1() {
  blaze::DynamicVector<int> dv_int(3);
  dv_int[0] = 1;
  dv_int[1] = 2;
  dv_int[2] = 4;

  return Rcpp::List::create(
    _["test"] = dv_int
  );
}


/*
Rcpp::List QRsolver(const blaze::DynamicMatrix<double>& X, const blaze::DynamicVector<double>& y) {
  blaze::DynamicMatrix<double> Q;
  blaze::DynamicMatrix<double> R;
  qr(X, Q, R);

  blaze::DynamicMatrix<double> S(R);
  blaze::invert(S);
  blaze::DynamicVector<double> coef(S * blaze::trans(Q) * y);

  blaze::DynamicVector<double> fitted = X * coef;
  blaze::DynamicVector<double> resid  = y - fitted;
  double s = std::sqrt((resid, resid ) / ((double) X.rows() - (double) R.rows()));

  blaze::DynamicVector<double> se(S.rows());
  for (size_t i=0UL; i<S.rows(); ++i) {
    se[i] = std::sqrt((row(S, i), row(S, i))) * s;
  }

  return List::create(
    _["coefficients"]  = coef,
    _["se"]            = se,
    _["rank"]          = (unsigned int) R.rows(),
    _["df.residual"]   = (unsigned int) X.rows() - (unsigned int) R.rows(),
    _["residuals"]     = resid,
    _["s"]             = s,
    _["fitted.values"] = fitted
  );
}

Rcpp::List LDLTSolver(const blaze::DynamicMatrix<double>& X, const blaze::DynamicVector<double>& y) {
  blaze::DynamicMatrix<double> XTXinv(blaze::trans(X) * X);
  blaze::invert<blaze::byLDLT>(XTXinv);

  blaze::DynamicVector<double> coef(XTXinv * blaze::trans(X) * y);

  blaze::DynamicVector<double> fitted = X * coef;
  blaze::DynamicVector<double> resid  = y - fitted;
  double s = std::sqrt((resid, resid) / ((double) X.rows() - (double) XTXinv.columns()));

  blaze::DynamicVector<double> se(XTXinv.columns());
  for (size_t i=0UL; i<XTXinv.columns(); ++i) {
    se[i] = std::sqrt( XTXinv(i, i) ) * s;
  }

  return List::create(
    _["coefficients"]  = coef,
    _["se"]            = se,
    _["rank"]          = (unsigned int) XTXinv.columns(),
    _["df.residual"]   = (unsigned int) X.rows() - (unsigned int) XTXinv.columns(),
    _["residuals"]     = resid,
    _["s"]             = s,
    _["fitted.values"] = fitted
  );
}

Rcpp::List LLTSolver(const blaze::DynamicMatrix<double>& X, const blaze::DynamicVector<double>& y) {
  blaze::DynamicMatrix<double> XTXinv(blaze::trans(X) * X);
  blaze::invert<blaze::byLLH>(XTXinv);

  blaze::DynamicVector<double> coef(XTXinv * blaze::trans(X) * y);

  blaze::DynamicVector<double> fitted = X * coef;
  blaze::DynamicVector<double> resid  = y - fitted;
  double s = std::sqrt( ( resid, resid ) / ((double) X.rows() - (double) XTXinv.columns()));

  blaze::DynamicVector<double> se(XTXinv.columns());
  for (size_t i=0UL; i<XTXinv.columns(); ++i) {
    se[i] = std::sqrt(XTXinv(i, i)) * s;
  }

  return List::create(
    _["coefficients"]  = coef,
    _["se"]            = se,
    _["rank"]          = (unsigned int) XTXinv.columns(),
    _["df.residual"]   = (unsigned int) X.rows() - (unsigned int) XTXinv.columns(),
    _["residuals"]     = resid,
    _["s"]             = s,
    _["fitted.values"] = fitted
  );
}

//' linear model fitting function based on RcppBlaze
//'
//' \code{fastLmPure} provides the estimates of the linear model based on \code{RcppBlaze}.
//'
//' \code{fastLm} estimates the linear model using the \code{solve}.
//'
//' @param X A model matrix.
//' @param y A response vector.
//' @param type A integer. 0 is QR solver, 1 is LLT solver and 2 is LDLT sovler.
//' @return A list containing coefficients, standard errors, rank of model matrix,
//'   degree of freedom of residuals, residuals, the standard deviation of random errors and
//'   fitted values.
//' @examples
//' # according to fastLm example in RcppArmadillo
//' data(trees, package="datasets")
//' flm <- fastLmPure(cbind(1, log(trees$Girth)), log(trees$Volume), 0)
//' print(flm)
//' @export
// [[Rcpp::export]]
List fastLmPure(blaze::DynamicMatrix<double> X, blaze::DynamicVector<double> y, int type) {
  if (X.rows() != y.size()) {
    throw std::invalid_argument("size mismatch");
  }

  switch(type) {
    case QRSolverType:
      return QRsolver(X, y);
    case LLTSolverType:
      return LLTSolver(X, y);
    case LDLTSolverType:
      return LDLTSolver(X, y);
    default:
      throw std::invalid_argument("invalid type");
  }
}

 */

