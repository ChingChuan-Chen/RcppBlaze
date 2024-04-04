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
Rcpp::List testAs1(Rcpp::List input_list) {
  blaze::DynamicMatrix<int, blaze::columnMajor> dm_int = input_list[0];
  Rcpp::Rcout << dm_int << std::endl;

  blaze::DynamicMatrix<double, blaze::columnMajor> dm_double = input_list[1];
  Rcpp::Rcout << dm_double  << std::endl;

  blaze::StaticMatrix<double, 2, 3, blaze::columnMajor> sm_double = input_list[1];
  Rcpp::Rcout << sm_double  << std::endl;

  return Rcpp::List::create(
    Rcpp::_["blaze::dm_int"] = blaze::sum(dm_int),
    Rcpp::_["blaze::dm_double"] = blaze::sum(dm_double),
    Rcpp::_["blaze::sm_double"] = blaze::sum(sm_double)
  );
}

// [[Rcpp::export]]
Rcpp::List testWrap1() {
  using namespace std::complex_literals;

  blaze::DynamicMatrix<int, blaze::columnMajor> dm_int_cm(2, 3);
  dm_int_cm = { { 3,  6, 3 }, { -4, -12, 0 } };

  blaze::DynamicMatrix<int, blaze::rowMajor> dm_int_rm(2, 3);
  dm_int_rm = { { 3,  6, 3 }, { -4, -12, 0 } };

  /*
  blaze::DynamicMatrix<std::complex<double>> dm_cmplx_cm(2, 3);
  dm_cmplx_cm = { {  1.5 + 0.2i,  2.5 + 0.4i, 4.5 - 0.5i }, { -0.5 - 0.1i, -12.1 + 2i, 0.6i } };

  blaze::DynamicMatrix<std::complex<double>, blaze::rowMajor> dm_cmplx_rm(2, 3);
  dm_cmplx_rm = { {  1.5 + 0.2i,  2.5 + 0.4i, 4.5 - 0.5i }, { -0.5 - 0.1i, -12.1 + 2i, 0.6i } };
  */

  blaze::DynamicMatrix<double, blaze::columnMajor> dm_double_cm(2, 3);
  dm_double_cm = { { 1.5,  2.5, 4.5 }, { -0.5, -12.1, -3.3 } };

  blaze::DynamicMatrix<double, blaze::rowMajor> dm_double_rm(2, 3);
  dm_double_rm = { { 1.5,  2.5, 4.5 }, { -0.5, -12.1, -3.3 } };

  blaze::StaticMatrix<double, 2UL, 3UL, blaze::columnMajor> sm_double_cm;
  sm_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::StaticMatrix<double, 2UL, 3UL, blaze::rowMajor> sm_double_rm;
  sm_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::HybridMatrix<double, 2UL, 3UL, blaze::columnMajor> hm_double_cm(2UL, 3UL);
  hm_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  blaze::HybridMatrix<double, 2UL, 3UL, blaze::rowMajor> hm_double_rm(2UL, 3UL);
  hm_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::Deallocate> al_pa_memory_cm(blaze::allocate<double>(16UL));
  blaze::CustomMatrix<double, blaze::aligned, blaze::padded, blaze::columnMajor> cm_al_pa_double_cm(al_pa_memory_cm.get(), 2UL, 3UL, 8UL);
  cm_al_pa_double_cm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  std::unique_ptr<double[], blaze::Deallocate> al_pa_memory_rm(blaze::allocate<double>(16UL));
  blaze::CustomMatrix<double, blaze::aligned, blaze::padded, blaze::rowMajor> cm_al_pa_double_rm(al_pa_memory_rm.get(), 2UL, 3UL, 8UL);
  cm_al_pa_double_rm = { { 1.5, -2.5, 4.5 }, { -5.5, 8.5, -7.3 } };

  return Rcpp::List::create(
    _["dm_int_cm"] = dm_int_cm,
    _["dm_int_rm"] = dm_int_rm,
   // _["dm_cmplx_cm"] = dm_cmplx_cm,
   // _["dm_cmplx_rm"] = dm_cmplx_rm,
    _["dm_double_cm"] = dm_double_cm,
    _["dm_double_rm"] = dm_double_rm,
    _["sm_double_cm"] = sm_double_cm,
    _["sm_double_rm"] = sm_double_rm,
    _["hm_double_cm"] = hm_double_cm,
    _["hm_double_rm"] = hm_double_rm,
    _["cm_al_pa_double_cm"] = cm_al_pa_double_cm,
    _["cm_al_pa_double_rm"] = cm_al_pa_double_rm
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

