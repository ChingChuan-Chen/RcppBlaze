// Copyright (C)  2017 - 2024  Ching-Chuan Chen
// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is based on fastLm.cpp and fastLm.h from RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.
//
// fastLm.cpp: lm model written RcppBlaze
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppBlaze is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppBlaze.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppBlaze.h>

using Rcpp::List;
using Rcpp::_;

enum {QRSolverType = 0, LDLTSolverType, LLTSolverType};

// [[Rcpp::export]]
Rcpp::List testAs1(Rcpp::NumericVector x) {
  blaze::DynamicVector<double, blaze::columnVector> y = Rcpp::as<blaze::DynamicVector<double, blaze::columnVector>>(x);
  Rcpp::Rcout << y << std::endl;
  return List::create(
    _["test"] = true
  );
}

// [[Rcpp::export]]
Rcpp::List testAs2(Rcpp::NumericVector x) {
  blaze::DynamicVector<double, blaze::rowVector> y = Rcpp::as<blaze::DynamicVector<double, blaze::rowVector>>(x);
  Rcpp::Rcout << y << std::endl;
  return List::create(
    _["test"] = true
  );
}

// [[Rcpp::export]]
Rcpp::List testAs3(Rcpp::NumericMatrix x) {
  blaze::DynamicMatrix<double, blaze::columnMajor> y = Rcpp::as<blaze::DynamicMatrix<double, blaze::columnMajor>>(x);
  Rcpp::Rcout << y << std::endl;
  return List::create(
    _["test"] = true
  );
}

// [[Rcpp::export]]
Rcpp::List testAs4(Rcpp::NumericMatrix x) {
  blaze::DynamicMatrix<double, blaze::rowMajor> y = Rcpp::as<blaze::DynamicMatrix<double, blaze::rowMajor>>(x);
  Rcpp::Rcout << y << std::endl;
  return List::create(
    _["test"] = true
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

