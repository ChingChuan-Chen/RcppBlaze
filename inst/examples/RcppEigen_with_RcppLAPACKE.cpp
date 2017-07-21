// Copyright (C)  2017         Chingchuan Chen
// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is based on fastLm.cpp and fastLm.h from RcppEigen.
// This file is part of RcppBlaze.
//
// RcppBlaze_with_RcppLAPACKE.cpp: lm model written RcppBlaze
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

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
// must include RcppLAPACKE first
#include <lapacke.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen, RcppLAPACKE)]]

namespace lmsol {
using Eigen::ArrayXd;
using Eigen::ColPivHouseholderQR;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
using Eigen::HouseholderQR;
using Eigen::JacobiSVD;
using Eigen::LDLT;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::SelfAdjointView;
using Eigen::TriangularView;
using Eigen::VectorXd;
using Eigen::Upper;

using Rcpp::_;
using Rcpp::as;
using Rcpp::CharacterVector;
using Rcpp::clone;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::RObject;
using Rcpp::wrap;

using std::invalid_argument;
using std::numeric_limits;

typedef MatrixXd::Index                                 Index;
typedef MatrixXd::Scalar                                Scalar;
typedef MatrixXd::RealScalar                            RealScalar;
typedef ColPivHouseholderQR<MatrixXd>::PermutationType  Permutation;
typedef Permutation::IndicesType                        Indices;

class lm {
protected:
  Map<MatrixXd> m_X;	/**< model matrix */
  Map<VectorXd> m_y;	/**< response vector */
  Index         m_n;	/**< number of rows of X */
  Index         m_p;	/**< number of columns of X */
  VectorXd      m_coef;	/**< coefficient vector */
  int           m_r;	/**< computed rank or NA_INTEGER */
  VectorXd      m_fitted;	/**< vector of fitted values */
  VectorXd      m_se;	/**< standard errors  */
  RealScalar    m_prescribedThreshold; /**< user specified tolerance */
  bool          m_usePrescribedThreshold;
public:
  lm(const Map<MatrixXd>&, const Map<VectorXd>&);

  ArrayXd                Dplus(const ArrayXd& D);
  MatrixXd                 I_p() const {return MatrixXd::Identity(m_p, m_p);}
  MatrixXd                 XtX() const;

  // setThreshold and threshold are based on ColPivHouseholderQR methods
  RealScalar         threshold() const;
  const VectorXd&           se() const {return m_se;}
  const VectorXd&         coef() const {return m_coef;}
  const VectorXd&       fitted() const {return m_fitted;}
  int                     rank() const {return m_r;}
  lm&             setThreshold(const RealScalar&);
};

class ColPivQR : public lm {
public:
  ColPivQR(const Map<MatrixXd>&, const Map<VectorXd>&);
};

class Llt : public lm {
public:
  Llt(const Map<MatrixXd>&, const Map<VectorXd>&);
};

class Ldlt : public lm {
public:
  Ldlt(const Map<MatrixXd>&, const Map<VectorXd>&);
};

class QR : public lm {
public:
  QR(const Map<MatrixXd>&, const Map<VectorXd>&);
};

class GESDD : public lm {
public:
  GESDD(const Map<MatrixXd>&, const Map<VectorXd>&);
};

class SVD : public lm {
public:
  SVD(const Map<MatrixXd>&, const Map<VectorXd>&);
};

class SymmEigen : public lm {
public:
  SymmEigen(const Map<MatrixXd>&, const Map<VectorXd>&);
};

lm::lm(const Map<MatrixXd> &X, const Map<VectorXd> &y)
  : m_X(X),
    m_y(y),
    m_n(X.rows()),
    m_p(X.cols()),
    m_coef(VectorXd::Constant(m_p, ::NA_REAL)),
    m_r(::NA_INTEGER),
    m_fitted(m_n),
    m_se(VectorXd::Constant(m_p, ::NA_REAL)),
    m_usePrescribedThreshold(false) {
}

lm& lm::setThreshold(const RealScalar& threshold) {
  m_usePrescribedThreshold = true;
  m_prescribedThreshold = threshold;
  return *this;
}

inline ArrayXd lm::Dplus(const ArrayXd& d) {
  ArrayXd   di(d.size());
  double  comp(d.maxCoeff() * threshold());
  for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1./d[j];
  m_r          = (di != 0.).count();
  return di;
}

MatrixXd lm::XtX() const {
  return MatrixXd(m_p, m_p).setZero().selfadjointView<Lower>().
  rankUpdate(m_X.adjoint());
}

/** Returns the threshold that will be used by certain methods such as rank().
*
*  The default value comes from experimenting (see "LU precision
*  tuning" thread on the Eigen list) and turns out to be
*  identical to Higham's formula used already in LDLt.
*
*  @return The user-prescribed threshold or the default.
*/
RealScalar lm::threshold() const {
  return m_usePrescribedThreshold ? m_prescribedThreshold
  : numeric_limits<double>::epsilon() * m_p;
}

ColPivQR::ColPivQR(const Map<MatrixXd> &X, const Map<VectorXd> &y)
  : lm(X, y) {
  ColPivHouseholderQR<MatrixXd> PQR(X); // decompose the model matrix
  Permutation                  Pmat(PQR.colsPermutation());
  m_r                               = PQR.rank();
  if (m_r == m_p) {	// full rank case
    m_coef     = PQR.solve(y);
    m_fitted   = X * m_coef;
    m_se       = Pmat * PQR.matrixQR().topRows(m_p).
    triangularView<Upper>().solve(I_p()).rowwise().norm();
    return;
  }
  MatrixXd                     Rinv(PQR.matrixQR().topLeftCorner(m_r, m_r).
                                      triangularView<Upper>().
                                      solve(MatrixXd::Identity(m_r, m_r)));
  VectorXd                  effects(PQR.householderQ().adjoint() * y);
  m_coef.head(m_r)                  = Rinv * effects.head(m_r);
  m_coef                            = Pmat * m_coef;
  // create fitted values from effects
  // (can't use X*m_coef if X is rank-deficient)
  effects.tail(m_n - m_r).setZero();
  m_fitted                          = PQR.householderQ() * effects;
  m_se.head(m_r)                    = Rinv.rowwise().norm();
  m_se                              = Pmat * m_se;
}

QR::QR(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
  HouseholderQR<MatrixXd> QR(X);
  m_coef                     = QR.solve(y);
  m_fitted                   = X * m_coef;
  m_se                       = QR.matrixQR().topRows(m_p).
  triangularView<Upper>().solve(I_p()).rowwise().norm();
}


Llt::Llt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
  LLT<MatrixXd>  Ch(XtX().selfadjointView<Lower>());
  m_coef            = Ch.solve(X.adjoint() * y);
  m_fitted          = X * m_coef;
  m_se              = Ch.matrixL().solve(I_p()).colwise().norm();
}

Ldlt::Ldlt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
  LDLT<MatrixXd> Ch(XtX().selfadjointView<Lower>());
  Dplus(Ch.vectorD());	// to set the rank
  //FIXME: Check on the permutation in the LDLT and incorporate it in
  //the coefficients and the standard error computation.
  //	m_coef            = Ch.matrixL().adjoint().
  //	    solve(Dplus(D) * Ch.matrixL().solve(X.adjoint() * y));
  m_coef            = Ch.solve(X.adjoint() * y);
  m_fitted          = X * m_coef;
  m_se              = Ch.solve(I_p()).diagonal().array().sqrt();
}

int gesdd(MatrixXd& A, ArrayXd& S, MatrixXd& Vt) {
  int info, mone = -1, m = A.rows(), n = A.cols();
  std::vector<int> iwork(8 * n);
  double wrk;
  if (m < n || S.size() != n || Vt.rows() != n || Vt.cols() != n)
    throw std::invalid_argument("dimension mismatch in gesvd");
  F77_CALL(dgesdd)((char*) "O", &m, &n, A.data(), &m, S.data(), A.data(),
           &m, Vt.data(), &n, &wrk, &mone, &iwork[0], &info);
  int lwork(wrk);
  std::vector<double> work(lwork);
  F77_CALL(dgesdd)((char*) "O", &m, &n, A.data(), &m, S.data(), A.data(),
           &m, Vt.data(), &n, &work[0], &lwork, &iwork[0], &info);
  return info;
}

GESDD::GESDD(const Map<MatrixXd>& X, const Map<VectorXd> &y) : lm(X, y) {
  MatrixXd   U(X), Vt(m_p, m_p);
  ArrayXd   S(m_p);
  if (gesdd(U, S, Vt)) throw std::runtime_error("error in gesdd");
  MatrixXd VDi(Vt.adjoint() * Dplus(S).matrix().asDiagonal());
  m_coef      = VDi * U.adjoint() * y;
  m_fitted    = X * m_coef;
  m_se        = VDi.rowwise().norm();
}

SVD::SVD(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
  JacobiSVD<MatrixXd>  UDV(X.jacobiSvd(ComputeThinU|ComputeThinV));
  MatrixXd             VDi(UDV.matrixV() *
    Dplus(UDV.singularValues().array()).matrix().asDiagonal());
  m_coef                   = VDi * UDV.matrixU().adjoint() * y;
  m_fitted                 = X * m_coef;
  m_se                     = VDi.rowwise().norm();
}

SymmEigen::SymmEigen(const Map<MatrixXd> &X, const Map<VectorXd> &y)
  : lm(X, y) {
  SelfAdjointEigenSolver<MatrixXd> eig(XtX().selfadjointView<Lower>());
  MatrixXd   VDi(eig.eigenvectors() *
    Dplus(eig.eigenvalues().array()).sqrt().matrix().asDiagonal());
  m_coef         = VDi * VDi.adjoint() * X.adjoint() * y;
  m_fitted       = X * m_coef;
  m_se           = VDi.rowwise().norm();
}

enum {ColPivQR_t = 0, QR_t, LLT_t, LDLT_t, SVD_t, SymmEigen_t, GESDD_t};

static inline lm do_lm(const Map<MatrixXd> &X, const Map<VectorXd> &y, int type) {
  switch(type) {
  case ColPivQR_t:
    return ColPivQR(X, y);
  case QR_t:
    return QR(X, y);
  case LLT_t:
    return Llt(X, y);
  case LDLT_t:
    return Ldlt(X, y);
  case SVD_t:
    return SVD(X, y);
  case SymmEigen_t:
    return SymmEigen(X, y);
  case GESDD_t:
    return GESDD(X, y);
  }
  throw invalid_argument("invalid type");
  return ColPivQR(X, y);	// -Wall
}

List fastLm(Rcpp::NumericMatrix Xs, Rcpp::NumericVector ys, int type) {
  const Map<MatrixXd>  X(as<Map<MatrixXd> >(Xs));
  const Map<VectorXd>  y(as<Map<VectorXd> >(ys));
  Index                n = X.rows();
  if ((Index)y.size() != n) throw invalid_argument("size mismatch");

  // Select and apply the least squares method
  lm                 ans(do_lm(X, y, type));

  // Copy coefficients and install names, if any
  NumericVector     coef(wrap(ans.coef()));

  List          dimnames(NumericMatrix(Xs).attr("dimnames"));
  if (dimnames.size() > 1) {
    RObject   colnames = dimnames[1];
    if (!(colnames).isNULL())
      coef.attr("names") = clone(CharacterVector(colnames));
  }

  VectorXd         resid = y - ans.fitted();
  int               rank = ans.rank();
  int                 df = (rank == ::NA_INTEGER) ? n - X.cols() : n - rank;
  double               s = resid.norm() / std::sqrt(double(df));
  // Create the standard errors
  VectorXd            se = s * ans.se();

  return List::create(_["coefficients"]  = coef,
                      _["se"]            = se,
                      _["rank"]          = rank,
                      _["df.residual"]   = df,
                      _["residuals"]     = resid,
                      _["s"]             = s,
                      _["fitted.values"] = ans.fitted());
}
}

// This defines the R-callable function 'fastLm'
// [[Rcpp::export]]
Rcpp::List fastLm_Impl2(Rcpp::NumericMatrix X, Rcpp::NumericVector y, int type) {
  return lmsol::fastLm(X, y, type);
}
