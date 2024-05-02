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
#if BLAZE_OPENMP_PARALLEL_MODE
#include <omp.h>
#endif

using Rcpp::_;

typedef typename blaze::CustomVector<double, blaze::aligned, blaze::padded> BlazeDblCuVector;
typedef typename blaze::CustomMatrix<double, blaze::aligned, blaze::padded, blaze::columnMajor> BlazeDblCuMatrix;

#define INIT_VEC(__var__, __data__, __size__, __padded_size__)                                     \
std::unique_ptr<double[], blaze::Deallocate> __data__( blaze::allocate<double>(__padded_size__) ); \
BlazeDblCuVector __var__ ( __data__.get(), __size__, __padded_size__ );\

#define INIT_MAT(__var__, __data__, __rows__, __cols__, __padded__)                                      \
std::unique_ptr<double[], blaze::Deallocate> __data__( blaze::allocate<double>(__padded__ * __cols__) ); \
BlazeDblCuMatrix __var__ ( __data__.get(), __rows__, __cols__, __padded__ );

enum {QRSolverType = 0, LDLTSolverType, LLTSolverType};

Rcpp::List QRsolver(const BlazeDblCuMatrix& X, const BlazeDblCuVector& y, size_t n_padded, size_t p_padded) {
  const size_t n = X.rows(), p = X.columns();

  INIT_MAT(Q, q_data, n, p, n_padded);
  INIT_MAT(R, r_data, p, p, p_padded);
  blaze::qr(X, Q, R);

  const size_t rank = R.rows();
  INIT_VEC(coef, coef_data, p, p_padded);
  blaze::invert(R);
  coef = R * blaze::trans(Q) * y;

  INIT_VEC(fitted, fitted_data, n, n_padded);
  fitted = X * coef;
  INIT_VEC(resid, resid_data, n, n_padded);
  resid = y - fitted;
  double s = std::sqrt(blaze::dot(resid, resid) / ((double) (n-p)));
  INIT_VEC(se, se_data, p, p_padded);
  for (size_t i=0UL; i<p; ++i) {
    se[i] = std::sqrt(blaze::dot(row(R, i), row(R, i))) * s;
  }

  return Rcpp::List::create(
    _["coefficients"]  = coef,
    _["se"]            = se,
    _["rank"]          = (int)rank,
    _["df.residual"]   = (int)(n-rank),
    _["residuals"]     = resid,
    _["s"]             = s,
    _["fitted.values"] = fitted
  );
}

Rcpp::List LLTSolver(const BlazeDblCuMatrix& X, const BlazeDblCuVector& y, size_t n_padded, size_t p_padded) {
  const size_t n = X.rows(), p = X.columns();

  INIT_MAT(XTXinv, xtx_data, p, p, p_padded);
  XTXinv = blaze::trans(X) * X;
  blaze::invert<blaze::byLLH>(XTXinv);

  INIT_VEC(coef, coef_data, p, p_padded);
  coef = XTXinv * blaze::trans(X) * y;

  INIT_VEC(fitted, fitted_data, n, n_padded);
  fitted = X * coef;
  INIT_VEC(resid, resid_data, n, n_padded);
  resid = y - fitted;
  double s = std::sqrt(blaze::dot(resid, resid) / ((double) (n-p)));
  INIT_VEC(se, se_data, p, p_padded);
  se = blaze::diagonal(XTXinv) * s;

  return Rcpp::List::create(
    _["coefficients"]  = coef,
    _["se"]            = se,
    _["rank"]          = (int) p,
    _["df.residual"]   = (int) (n-p),
    _["residuals"]     = resid,
    _["s"]             = s,
    _["fitted.values"] = fitted
  );
}

Rcpp::List LDLTSolver(const BlazeDblCuMatrix& X, const BlazeDblCuVector& y, size_t n_padded, size_t p_padded) {
  const size_t n = X.rows(), p = X.columns();

  INIT_MAT(XTXinv, xtx_data, p, p, p_padded);
  XTXinv = blaze::trans(X) * X;
  blaze::invert<blaze::byLDLT>(XTXinv);

  INIT_VEC(coef, coef_data, p, p_padded);
  coef = XTXinv * blaze::trans(X) * y;

  INIT_VEC(fitted, fitted_data, n, n_padded);
  fitted = X * coef;
  INIT_VEC(resid, resid_data, n, n_padded);
  resid = y - fitted;
  double s = std::sqrt(blaze::dot(resid, resid) / ((double) (n-p)));
  INIT_VEC(se, se_data, p, p_padded);
  se = blaze::diagonal(XTXinv) * s;

  return Rcpp::List::create(
    _["coefficients"]  = coef,
    _["se"]            = se,
    _["rank"]          = (int) p,
    _["df.residual"]   = (int) (n-p),
    _["residuals"]     = resid,
    _["s"]             = s,
    _["fitted.values"] = fitted
  );
}

//' linear model fitting function based on RcppBlaze
//'
//' \code{fastLmPure} provides the estimates of the linear model based on \strong{RcppBlaze}.
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
Rcpp::List fastLmPure(Rcpp::NumericMatrix X, Rcpp::NumericVector y, int type) {
  if (X.nrow() != y.size()) {
    throw std::invalid_argument("size mismatch");
  }

  // define sizes
  const size_t n = (size_t) X.nrow(), p = (size_t) X.ncol();
  const std::size_t n_padded = blaze::nextMultiple<std::size_t>(n, blaze::SIMDTrait<double>::size);
  const std::size_t p_padded = blaze::nextMultiple<std::size_t>(p, blaze::SIMDTrait<double>::size);

  // define y
  INIT_VEC(y_, y_data, n, n_padded);
  RcppBlaze::copyToCustomVector(y, y_);

  // define x
  INIT_MAT(x_, x_data, n, p, n_padded);
  RcppBlaze::copyToCustomMatrix(X, x_);

  switch(type) {
    case QRSolverType:
      return QRsolver(x_, y_, n_padded, p_padded);
    case LLTSolverType:
      return LLTSolver(x_, y_, n_padded, p_padded);
    case LDLTSolverType:
      return LDLTSolver(x_, y_, n_padded, p_padded);
    default:
      throw std::invalid_argument("invalid type");
  }
}
