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
Rcpp::List testAs1(Rcpp::List input_list) {
  blaze::CompressedVector<double, blaze::rowVector> y = input_list[0];
  Rcpp::Rcout << "CompressedVector rowVector:" << std::endl << y << std::endl;

  blaze::CompressedVector<double, blaze::columnVector> x = input_list[1];
  Rcpp::Rcout << "CompressedVector columnVector:" << std::endl << x << std::endl;

  blaze::ZeroVector<double> zv_double_dgCMatrix = input_list[2];
  Rcpp::Rcout << "ZeroVector columnVector:" << std::endl << zv_double_dgCMatrix << std::endl;

  return Rcpp::List::create(
    _["test"] = true
  );
}


using namespace std::complex_literals;

// [[Rcpp::export]]
Rcpp::List testWrap1() {
  blaze::UniformMatrix<double, blaze::columnMajor> um_double_cm(2UL, 3UL, -3.2);
  blaze::UniformMatrix<double, blaze::rowMajor> um_double_rm(2UL, 3UL, -3.2);
  blaze::UniformMatrix<std::complex<double>, blaze::columnMajor> um_cplx(2UL, 3UL, -1.8 + 0.6i);

  blaze::IdentityMatrix<double, blaze::columnMajor> im_double(3UL);

  blaze::ZeroMatrix<double, blaze::columnMajor> zm_double_cm(2UL, 3UL);
  blaze::ZeroMatrix<double, blaze::rowMajor> zm_double_rm(2UL, 3UL);
  return Rcpp::List::create(
    _["um_double_cm"] = um_double_cm,
    _["um_double_rm"] = um_double_rm,
    _["um_cplx"] = um_cplx,
    _["im_double"] = im_double,
    _["zm_double_cm"] = zm_double_cm,
    _["zm_double_rm"] = zm_double_rm
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

