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

// [[Rcpp::export]]
Rcpp::List testWrap() {
  blaze::StaticMatrix<double, 3UL, 3UL> D{
    {  3.0,  0.0,  0.0 },
    {  8.0,  0.0,  0.0 },
    { -2.0, -1.0,  4.0 }
  };

  auto rs1 = blaze::rows<0UL, 2UL>(D);
  auto rs2 = blaze::rows(D, {0UL, 2UL});

  blaze::CompressedMatrix<double> E(3UL, 3UL);
  E.reserve(5);
  E.append(0, 0, 3.0);
  E.finalize(0);
  E.append(1, 0, 8.0);
  E.finalize(1);
  E.append(2, 0, -2.0);
  E.append(2, 1, -1.0);
  E.append(2, 2, 4.0);
  E.finalize(2);

  auto rs3 = blaze::rows<0UL, 2UL>(E);
  auto rs4 = blaze::rows(E, {0UL, 2UL});

  return Rcpp::List::create(
    _["row0"] = rs1,
    _["row1"] = rs2,
    _["row2"] = rs3,
    _["row3"] = rs4
  );
/*
  blaze::StaticVector<double, 3UL> a;
  a = {1.5, -2.5, 4.5};

  blaze::StaticVector<double, 3UL> b;
  b = {-1.5, 2.5, 0.0};

  return Rcpp::List::create(
    _["2.0 + a"] = Rcpp::RcppBlaze::blaze_dv_expr_wrap(2.0 + a)
    _["2.0 - a"] = 2.0 - a,
    _["2.0 * a"] = 2.0 * a,
    _["a+b"] = a + b,
    _["a-b"] = a - b,
    _["a*b"] = a * b
  );
  */
}


enum {QRSolverType = 0, LDLTSolverType, LLTSolverType};

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

  return Rcpp::List::create(
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

  return Rcpp::List::create(
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

  return Rcpp::List::create(
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
Rcpp::List fastLmPure(blaze::DynamicMatrix<double> X, blaze::DynamicVector<double> y, int type) {
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
