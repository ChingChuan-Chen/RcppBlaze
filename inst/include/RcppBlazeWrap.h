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

namespace Rcpp {

  namespace RcppBlaze {

    template <typename VT, bool TF>
    SEXP blaze_wrap(const blaze::DenseVector<VT, TF>& x) {
      typedef typename VT::ElementType Type;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
	    size_t size = (*x).size();
	    Rcpp::Vector<RTYPE> out(size);
      for (size_t i=0UL; i<size; ++i) {
        out[i] = Rcpp::internal::caster<Type, value_t>((*x)[i]);
      }
      return out;
    }

    template <typename Type, bool TF>
    SEXP blaze_uv_wrap(const blaze::UniformVector<Type, TF>& x) {
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
	    size_t size = (*x).size();
      Rcpp::Vector<RTYPE> out(size);
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
          p[element->index() + 1UL] = 1;
          ++k;
        }
        for( size_t i=1UL; i<(size_t)p.size(); ++i) {
          p[i] = p[i-1UL] + p[i];
        }
      } else {
        p[1UL] = nonZeroSize;
        for (auto element=(*sv).begin(); element!=(*sv).end(); ++element) {
          x[k] = Rcpp::internal::caster<Type, double>(element->value());
          idx[k] = element->index();
          ++k;
        }
      }

      std::string klass = "dgCMatrix";
      Rcpp::S4 out(klass);
      if (TF == blaze::rowVector) {
        out.slot("Dim") = Rcpp::Dimension(1, (int)size);
      } else {
        out.slot("Dim") = Rcpp::Dimension((int)size, 1);
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
    SEXP blaze_wrap(const blaze::SparseMatrix<MT, SO>& sm) {
      typedef typename MT::ElementType Type;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      if ((RTYPE != INTSXP) && (RTYPE != REALSXP)) {
        Rcpp::stop("SparseMatrix only supports int, float and double!");
      }

      size_t nonZeroSize = (*sm).nonZeros(), m = (*sm).rows(), n = (*sm).columns();
      Rcpp::Vector<REALSXP> x(nonZeroSize);
      Rcpp::IntegerVector idx(nonZeroSize);
      size_t pSizeMinus1 = (SO == blaze::rowMajor)?m:n;
      Rcpp::IntegerVector p(pSizeMinus1 + 1UL);
      p[pSizeMinus1] = (int) nonZeroSize;

      size_t k = 0UL;
      for (size_t i=0UL; i<pSizeMinus1; ++i) {
        p[i] = (int) k;
        for (auto element=(*sm).begin(i); element!=(*sm).end(i); ++element) {
          x[k] = Rcpp::internal::caster<Type, double>(element->value());
          idx[k] = element->index();
          ++k;
        }
      }

      std::string klass = (SO == blaze::rowMajor)?"dgRMatrix":"dgCMatrix";
      Rcpp::S4 out(klass);
      out.slot("Dim") = Rcpp::Dimension(m, n);
      out.slot((SO == blaze::rowMajor)?"j":"i") = idx;
      out.slot("p") = p;
      out.slot("x") = x;
      return out;
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

  // wrap for adaptors for blaze matrices
  template <typename MT, bool SO, bool DF, bool NF>
  SEXP wrap(const blaze::SymmetricMatrix<MT, SO, DF, NF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

#define RCPPBLAZE_ADAPTOR_WRAPPER(__ADAPTOR__)                 \
  template <typename MT, bool SO, bool DF> SEXP wrap(          \
      const __ADAPTOR__<MT, SO, DF>& x                         \
  ) {return RcppBlaze::blaze_wrap(x);};                        \

RCPPBLAZE_ADAPTOR_WRAPPER(blaze::HermitianMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::LowerMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::UniLowerMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::StrictlyLowerMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::UpperMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::UniUpperMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::StrictlyUpperMatrix);
RCPPBLAZE_ADAPTOR_WRAPPER(blaze::DiagonalMatrix);

#undef RCPPBLAZE_ADAPTOR_WRAPPER

  // wrap for views for blaze vectors
  template <typename VT, blaze::AlignmentFlag AF, bool TF, bool DF, size_t... CSAs>
  SEXP wrap(const blaze::Subvector<VT, AF, TF, DF, CSAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, bool TF, bool DF, typename... CEAs>
  SEXP wrap(const blaze::Elements<VT, TF, DF, CEAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  // wrap for views for blaze matrices
  template <typename MT, bool SO, bool DF, bool SF, size_t... CRAs>
  SEXP wrap(const blaze::Row<MT, SO, DF, SF, CRAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, bool SO, bool DF, bool SF, size_t... CCAs>
  SEXP wrap(const blaze::Column<MT, SO, DF, SF, CCAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, bool TF, bool DF, bool MF, ptrdiff_t... CBAs>
  SEXP wrap(const blaze::Band<MT, TF, DF, MF, CBAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, blaze::AlignmentFlag AF, bool SO, bool DF, size_t... CSAs>
  SEXP wrap(const blaze::Submatrix<MT, AF, SO, DF, CSAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, bool SO, bool DF, bool SF, typename... CRAs>
  SEXP wrap(const blaze::Rows<MT, SO, DF, SF, CRAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, bool SO, bool DF, bool SF, typename... CCAs>
  SEXP wrap(const blaze::Columns<MT, SO, DF, SF, CCAs...>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  /* support blaze Expression */
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatAddExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatKronExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, typename OP, bool SO> SEXP wrap(const blaze::DMatDMatMapExpr<MT1, MT2, OP, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::DMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSchurExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSolveExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSubExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::DMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT, bool TF> SEXP wrap(const blaze::DMatDVecSolveExpr<MT, VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename OP, bool SO> SEXP wrap(const blaze::DMatMapExpr<MT, OP, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::DMatScalarDivExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::DMatScalarMultExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatAddExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatKronExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::DMatSMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::DMatSMatSchurExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatSubExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::DMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecAddExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecCrossExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecDivExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecKronExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, typename OP, bool TF> SEXP wrap(const blaze::DVecDVecMapExpr<VT1, VT2, OP, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, typename OP> SEXP wrap(const blaze::DVecDVecOuterExpr<VT1, VT2, OP>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename OP, bool TF> SEXP wrap(const blaze::DVecMapExpr<VT, OP, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::DVecScalarDivExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::DVecScalarMultExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecAddExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecCrossExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecKronExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2> SEXP wrap(const blaze::DVecSVecOuterExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::SMatDMatKronExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::SMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatDMatSchurExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::SMatDMatSubExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::SMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename OP, bool SO> SEXP wrap(const blaze::SMatMapExpr<MT, OP, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::SMatScalarDivExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::SMatScalarMultExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatKronExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatSchurExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::SMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecCrossExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecDivExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecKronExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2> SEXP wrap(const blaze::SVecDVecOuterExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename OP, bool TF> SEXP wrap(const blaze::SVecMapExpr<VT, OP, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::SVecScalarDivExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::SVecScalarMultExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecAddExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecCrossExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecKronExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2> SEXP wrap(const blaze::SVecSVecOuterExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TDMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::TDMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TDMatSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TDMatSMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TDMatSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::TDMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename MT> SEXP wrap(const blaze::TDVecDMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename MT> SEXP wrap(const blaze::TDVecSMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TSMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatDMatSchurExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatDMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::TSMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatKronExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatSchurExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename MT, typename VT> SEXP wrap(const blaze::TSMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename MT> SEXP wrap(const blaze::TSVecDMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template <typename VT, typename MT> SEXP wrap(const blaze::TSVecSMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };
}

#endif
