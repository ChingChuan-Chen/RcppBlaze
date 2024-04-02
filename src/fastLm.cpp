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
using Rcpp::_;

enum {QRSolverType = 0, LDLTSolverType, LLTSolverType};

// [[Rcpp::export]]
Rcpp::List testAs1(Rcpp::NumericVector x) {
  blaze::DynamicVector<double, blaze::columnVector> y = Rcpp::as<blaze::DynamicVector<double, blaze::columnVector>>(x);
  Rcpp::Rcout << y << std::endl;
  return Rcpp::List::create(
    _["test"] = true
  );
}

// [[Rcpp::export]]
Rcpp::List testWrap1() {
  using namespace std::complex_literals;

  blaze::DynamicVector<int> dv_int(3);
  dv_int[0] = 1;
  dv_int[1] = 2;
  dv_int[2] = 4;

  blaze::DynamicVector<double> dv_double(3);
  dv_double[0] = 1.5;
  dv_double[1] = 2.5;
  dv_double[2] = 4.5;

  blaze::DynamicVector<std::complex<double>> dv_cmpl(3);
  dv_cmpl[0] = 1.5 + 0.2i;
  dv_cmpl[1] = 0.75 + 0.3i;
  dv_cmpl[2] = 3.0 + 0.1i;

  blaze::StaticVector<double, 3> sv_double;
  sv_double[0] = 1.5;
  sv_double[1] = 2.5;
  sv_double[2] = 4.5;

  blaze::HybridVector<double, 3> hv_double(3);
  hv_double[0] = 1.5;
  hv_double[1] = 2.5;
  hv_double[2] = 4.5;

  std::vector<double> vec(3UL);
  blaze::CustomVector<double, blaze::unaligned, blaze::unpadded> cv_ua_up_double(&vec[0], 3UL);
  cv_ua_up_double[0] = 1.5;
  cv_ua_up_double[1] = 2.5;
  cv_ua_up_double[2] = 4.5;

  std::unique_ptr<double[]> cv_ua_pa_double_mem(new double[4]);
  blaze::CustomVector<double, blaze::unaligned, blaze::padded> cv_ua_pa_double(cv_ua_pa_double_mem.get(), 3UL, 4UL);
  cv_ua_pa_double[0] = 1.5;
  cv_ua_pa_double[1] = 2.5;
  cv_ua_pa_double[2] = 4.5;


  std::unique_ptr<double[], blaze::Deallocate> cv_al_up_double_mem(blaze::allocate<double>(3UL));
  blaze::CustomVector<double, blaze::aligned, blaze::unpadded> cv_al_up_double(cv_al_up_double_mem.get(), 3UL);
  cv_al_up_double[0] = 1.5;
  cv_al_up_double[1] = 2.5;
  cv_al_up_double[2] = 4.5;

  std::unique_ptr<double[], blaze::Deallocate> cv_al_pa_double_mem(blaze::allocate<double>(4UL));
  blaze::CustomVector<double, blaze::aligned, blaze::padded> cv_al_pa_double(cv_al_pa_double_mem.get(), 3UL, 4UL);
  cv_al_pa_double[0] = 1.5;
  cv_al_pa_double[1] = 2.5;
  cv_al_pa_double[2] = 4.5;

  return Rcpp::List::create(
    _["dv_int"] = dv_int,
    _["dv_double"] = dv_double,
    _["dv_cmpl"] = dv_cmpl,
    _["sv_double"] = sv_double,
    _["hv_double"] = hv_double,
    _["cv_ua_up_double"] = cv_ua_up_double,
    _["cv_ua_pa_double"] = cv_ua_pa_double,
    _["cv_al_up_double"] = cv_al_up_double,
    _["cv_al_pa_double"] = cv_al_pa_double
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

