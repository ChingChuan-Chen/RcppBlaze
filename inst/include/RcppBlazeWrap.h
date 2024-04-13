// Copyright (C)  2017 - 2024  Ching-Chuan Chen
//
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.

#ifndef RcppBlaze__RcppBlazeWrap__h
#define RcppBlaze__RcppBlazeWrap__h

#include "RcppBlazeForward.h"
#include <Rcpp.h>

// most code are modified from the definition of matrix in blaze

namespace Rcpp {

  namespace RcppBlaze {

    template <typename VT, bool TF>
    SEXP blaze_wrap(const blaze::DenseVector<VT, TF>& x) {
      typedef typename VT::ElementType Type;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      const Type* data_pointer = blaze::data(*x);
      Rcpp::Vector<RTYPE> out = Rcpp::wrap(data_pointer, data_pointer + (*x).size());
      return out;
    }

    template <typename VT, bool TF>
    SEXP blaze_dv_expr_wrap(const blaze::Expression<blaze::DenseVector<VT, TF>>& x) {
      return blaze_wrap<VT, TF>(x);
    }

    template <typename Type, bool TF>
    SEXP blaze_uv_wrap(const blaze::UniformVector<Type, TF>& x) {
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
      Rcpp::Vector<RTYPE> out((*x).size());
      std::fill(out.begin(), out.end(), Rcpp::internal::caster<Type, value_t>((*x)[0UL]));
      return out;
    }

    template <typename MT, bool SO>
    SEXP blaze_wrap(const blaze::DenseMatrix<MT, SO>& x) {
      typedef typename MT::ElementType Type;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;

      size_t m = (*x).rows(), n = (*x).columns();
      Rcpp::Matrix<RTYPE> out(m, n);
      if (SO == blaze::rowMajor) {
        for (size_t i=0UL; i<m; ++i) {
          for (size_t j=0UL; j<n; ++j) {
            out(i, j) = Rcpp::internal::caster<Type, value_t>((*x)(i, j));
          }
        }
      } else {
        for (size_t j=0UL; j<n; ++j) {
          for (size_t i=0UL; i<m; ++i) {
            out(i, j) = Rcpp::internal::caster<Type, value_t>((*x)(i, j));
          }
        }
      }
      return out;
    }

    template <typename MT, bool SO>
    SEXP blaze_dm_expr_wrap(const blaze::Expression<blaze::DenseMatrix<MT, SO>>& x) {
      return blaze_wrap<MT, SO>(x);
    }

    template <typename Type, bool SO>
    SEXP blaze_um_wrap(const blaze::UniformMatrix<Type, SO>& x) {
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
      Rcpp::Matrix<RTYPE> out((*x).rows(), (*x).columns());
      std::fill(out.begin(), out.end(), Rcpp::internal::caster<Type, value_t>((*x)(0UL, 0UL)));
      return out;
    }

    template <typename VT, bool TF>
    SEXP blaze_wrap(const blaze::SparseVector<VT, TF>& sv) {
      typedef typename VT::ElementType Type;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      if ((RTYPE != INTSXP) && (RTYPE != REALSXP)) {
        Rcpp::stop("SparseVector only supports int, float and double!");
      }

      size_t nonZeroSize = (*sv).nonZeros(), size = (*sv).size();
      Rcpp::Vector<REALSXP> x(nonZeroSize);
      size_t pSize = 2UL;
      if (TF == blaze::rowVector) {
        pSize = size + 1UL;
      }
      Rcpp::IntegerVector p(pSize, 0);
      Rcpp::IntegerVector idx(nonZeroSize);

      size_t k = 0UL;
      if (TF == blaze::rowVector) {
        for (auto element=(*sv).begin(); element!=(*sv).end(); ++element) {
          x[k] = Rcpp::internal::caster<Type, double>(element->value());
          p[element->index() + 1] = 1;
          ++k;
        }
        for( size_t i=1UL; i<(size_t)p.size(); ++i) {
          p[i] = p[i-1UL] + p[i];
        }
      } else {
        p[1] = nonZeroSize;
        for (auto element=(*sv).begin(); element!=(*sv).end(); ++element) {
          x[k] = Rcpp::internal::caster<Type, double>(element->value());
          idx[k] = element->index();
          ++k;
        }
      }

      std::string klass = "dgCMatrix";
      Rcpp::S4 out(klass);
      if (TF == blaze::rowVector) {
        out.slot("Dim") = Rcpp::Dimension(1, size);
      } else {
        out.slot("Dim") = Rcpp::Dimension(size, 1);
      }
      out.slot("i") = idx;
      out.slot("p") = p;
      out.slot("x") = x;
      return out;
    }

    template <typename VT, bool TF>
    SEXP blaze_sv_expr_wrap(const blaze::Expression<blaze::SparseVector<VT, TF>>& x) {
      return blaze_wrap<VT, TF>(x);
    }

    template <typename MT, bool SO>
    SEXP blaze_wrap(const blaze::SparseMatrix<MT,SO>& sm) {
      typedef typename MT::ElementType Type;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      if ((RTYPE != INTSXP) && (RTYPE != REALSXP)) {
        Rcpp::stop("SparseMatrix only supports int, float and double!");
      }

      size_t m = (*sm).rows(), n = (*sm).columns;

      std::string klass = "dgCMatrix";
      Rcpp::S4 s(klass);
      s.slot("Dim") = Rcpp::Dimension(m, n);
      return s;

/*
      const size_t index((SO == blaze::rowMajor)?A.rows():A.columns());

      Rcpp::Vector<RTYPE> x(A.nonZeros());
      Rcpp::IntegerVector idx(A.nonZeros());
      Rcpp::IntegerVector p(index + 1UL);

      size_t k = 0UL;
      for (size_t i=0UL; i<index; ++i) {
        p[i] = (int) k;
        for (ConstIterator element=A.begin(i); element!=A.end(i); ++element) {
          x[k] = element->value();
          idx[k] = element->index();
          ++k;
        }
      }
      p[index] = (int) A.nonZeros();

      Rcpp::S4 s(klass);
      s.slot("Dim") = ::Rcpp::Dimension(A.rows(), A.columns());
      s.slot((SO == blaze::rowMajor)?"j":"i") = idx;
      s.slot("p") = p;
      s.slot("x") = x;
*/
    }

    template <typename MT, bool SO>
    SEXP blaze_sm_expr_wrap(const blaze::Expression<blaze::SparseMatrix<MT, SO>>& x) {
      return blaze_wrap<MT, SO>(x);
    }

  } // namespace RcppBlaze

  // wrap for blaze dense vector
  template <typename Type, bool TF>
  SEXP wrap(const blaze::DynamicVector<Type, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::StaticVector<Type, N, TF, AF, PF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::HybridVector<Type, N, TF, AF, PF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool TF>
  SEXP wrap(const blaze::CustomVector<Type, AF, PF, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, bool TF>
  SEXP wrap(const blaze::UniformVector<Type, TF>& x) {
    return RcppBlaze::blaze_uv_wrap(x);
  };

  // wrap for blaze dense matrix
  template <typename Type, bool SO>
  SEXP wrap(const blaze::DynamicMatrix<Type, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::StaticMatrix<Type, M, N, SO, AF, PF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::HybridMatrix<Type, M, N, SO, AF, PF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool SO>
  SEXP wrap(const blaze::CustomMatrix<Type, AF, PF, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, bool SO>
  SEXP wrap(const blaze::UniformMatrix<Type, SO>& x) {
    return RcppBlaze::blaze_um_wrap(x);
  };

  // wrap for blaze sparse vector
  template <typename Type, bool TF>
  SEXP wrap(const blaze::CompressedVector<Type, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, bool TF>
  SEXP wrap(const blaze::ZeroVector<Type, TF>& x) {
    return RcppBlaze::blaze_sv_expr_wrap(x);
  };

  // wrap for blaze sparse matrix
  template <typename Type, bool SO>
  SEXP wrap(const blaze::CompressedMatrix<Type, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename Type, bool SO>
  SEXP wrap(const blaze::IdentityMatrix<Type, SO>& x) {
    const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
    if ((RTYPE != INTSXP) && (RTYPE != REALSXP)) {
      Rcpp::stop("IdentityMatrix only supports int, float and double!");
    }
    size_t m = (*x).rows(), n = (*x).columns();
    if (m != n) {
      Rcpp::stop("IdentityMatrix must have the same numbers of rows and columns!");
    }
    std::string klass = "ddiMatrix";
    Rcpp::S4 out(klass);
    out.slot("Dim") = Rcpp::Dimension(m, n);
    out.slot("diag") = "N";
    out.slot("x") = Rcpp::Vector<RTYPE>(3, 1.0);
    return out;
  };

  template <typename Type, bool SO>
  SEXP wrap(const blaze::ZeroMatrix<Type, SO>& x) {
    const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
    if ((RTYPE != INTSXP) && (RTYPE != REALSXP)) {
      Rcpp::stop("ZeroMatrix only supports int, float and double!");
    }
    size_t m = (*x).rows(), n = (*x).columns();
    std::string klass = "dgCMatrix";
    Rcpp::S4 out(klass);
    out.slot("Dim") = Rcpp::Dimension(m, n);
    out.slot("i") = Rcpp::Vector<INTSXP>(0);
    out.slot("p") = Rcpp::Vector<INTSXP>(n+1, 0.0);
    out.slot("x") = Rcpp::Vector<RTYPE>(0);
    return out;
  };


  /* adaptors
   // TODO:
   // Convert Symmetric Matrices, Hermitian Matrices, Triangular Matrices
   // Convert Submatrix, Subvector, Row, Column, Rows, Columns, Band

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::DiagonalMatrix<MT, SO, DF>& dm) {
    return RcppBlaze::blaze_wrap(dm);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::LowerMatrix<MT, SO, DF>& lm) {
    return RcppBlaze::blaze_wrap(lm);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::UpperMatrix<MT, SO, DF>& um) {
    return RcppBlaze::blaze_wrap(um);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::HermitianMatrix<MT, SO, DF>& hm) {
    return RcppBlaze::blaze_wrap(hm);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::StrictlyLowerMatrix<MT, SO, DF>& slm) {
    return RcppBlaze::blaze_wrap(slm);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::StrictlyUpperMatrix<MT, SO, DF>& sum) {
    return RcppBlaze::blaze_wrap(sum);
  };

  template <typename MT, bool SO, bool DF, bool NF>
  SEXP wrap(const blaze::SymmetricMatrix<MT, SO, DF, NF>& sm) {
    return RcppBlaze::blaze_wrap(sm);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::UniLowerMatrix<MT, SO, DF>& ulm) {
    return RcppBlaze::blaze_wrap(ulm);
  };

  template <typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::UniUpperMatrix<MT, SO, DF>& uum) {
    return RcppBlaze::blaze_wrap(uum);
  };
   */

  // TODO:
  // Convert Exprs
}

#endif
