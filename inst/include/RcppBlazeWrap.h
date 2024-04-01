// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
// Copyright (C)  2017 - 2024  Ching-Chuan Chen
//
// This file is based on RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.
//
// RcppBlazeWrap.h: Rcpp/Blaze glue
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

#ifndef RcppBlaze__RcppBlazeWrap__h
#define RcppBlaze__RcppBlazeWrap__h

#include <RcppBlazeForward.h>
#include <Rcpp.h>

// most code are modified from the definition of matrix in blaze

// TODO:
// Convert IdentityMatrix, ZeroMatrix
// Convert Symmetric Matrices, Hermitian Matrices, Triangular Matrices
// Convert Submatrix, Subvector, Row, Column, Rows, Columns, Band

namespace Rcpp {

  namespace RcppBlaze {

    template<typename MT, bool SO>
    SEXP blaze_wrap(const blaze::DenseMatrix<MT, SO>& dm) {
      typedef typename MT::ElementType ET;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<ET>::rtype;
      ::Rcpp::RObject x = ::Rcpp::wrap(blaze::data(dm), dm.data() + dm.size());
      x.attr("dim") = ::Rcpp::Dimension(dm.rows(), dm.columns());
      return x;
    }

    template<typename VT, bool TF>
    SEXP blaze_wrap(const blaze::DenseVector<VT,TF>& dv) {
      typedef typename VT::ElementType ET;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<ET>::rtype;
      ::Rcpp::RObject x = ::Rcpp::wrap(blaze::data(dv), dv.data() + dv.size());
      return x;
    }


    /*
     template<typename MT, bool SO>
     SEXP blaze_wrap(const blaze::SparseMatrix<MT,SO>& sm)
     {
     typedef typename MT::ResultType     RT;
     typedef typename MT::ReturnType     RN;
     typedef typename MT::CompositeType  CT;
     typedef typename blaze::If<blaze::IsExpression<RN>, const RT, CT>::Type  Tmp;
     typedef typename blaze::RemoveReference<Tmp>::Type::ConstIterator ConstIterator;

     Tmp A(~sm);  // Evaluation of the sparse matrix operand

     const int RTYPE = ::Rcpp::traits::r_sexptype_traits<typename MT::ElementType>::rtype;
     std::string klass = "dgCMatrix";

     const size_t index((SO == blaze::rowMajor)?A.rows():A.columns());

     ::Rcpp::Vector<RTYPE> x(A.nonZeros());
     ::Rcpp::IntegerVector idx(A.nonZeros());
     ::Rcpp::IntegerVector p(index + 1UL);

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

     ::Rcpp::S4 s(klass);
     s.slot("Dim") = ::Rcpp::Dimension(A.rows(), A.columns());
     s.slot((SO == blaze::rowMajor)?"j":"i") = idx;
     s.slot("p") = p;
     s.slot("x") = x;
     return s;
     }

     template<typename VT, bool TF>
     SEXP blaze_wrap(const blaze::SparseVector<VT, TF>& sv)
     {
     typedef typename VT::ElementType    ET;
     typedef typename VT::CompositeType  CT;
     typedef typename blaze::RemoveReference<CT>::Type::ConstIterator  ConstIterator;

     CT a(~sv);  // Evaluation of the sparse vector operand

     const int RTYPE = ::Rcpp::traits::r_sexptype_traits<typename VT::ElementType>::rtype;
     std::string klass = "dgCMatrix";

     ::Rcpp::Vector<RTYPE> x(a.nonZeros());
     ::Rcpp::IntegerVector idx(a.nonZeros());
     size_t pSize = 2UL;
     if (TF == blaze::rowVector) {
     pSize = a.size() + 1UL;
     }
     ::Rcpp::IntegerVector p( pSize, 0 );

     size_t k = 0UL;
     if (TF == blaze::columnVector) {
     p[1UL] = a.nonZeros();
     for (ConstIterator element=a.begin(); element!=a.end(); ++element) {
     x[k] = element->value();
     idx[k] = element->index();
     ++k;
     }
     } else if (TF == blaze::rowVector) {
     for (ConstIterator element=a.begin(); element!=a.end(); ++element) {
     x[k] = element->value();
     p[element->index() + 1] = 1;
     ++k;
     }
     for( size_t i=1UL; i<(size_t)p.size();++i)
     p[i] = p[i-1UL] + p[i];
     }

     ::Rcpp::S4 s(klass);
     if (TF == blaze::rowVector) {
     s.slot("Dim") = ::Rcpp::Dimension(1, a.size());
     } else {
     s.slot("Dim") = ::Rcpp::Dimension(a.size(), 1);
     }
     s.slot("i") = idx;
     s.slot("p") = p;
     s.slot("x") = x;
     return s;
     }
     */
  } // namespace RcppBlaze

/* wrap for the child classes of DenseMatrix, DenseVector, SparseMatrix and SparseVector */

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DenseMatrix<MT, SO>& dm) {
    return RcppBlaze::blaze_wrap(dm);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::DenseVector<VT, TF>& dv) {
    return RcppBlaze::blaze_wrap(dv);
  };
  /*
   template<typename VT, bool TF>
   SEXP wrap(const blaze::SparseVector<VT, TF>& sv) {
   return RcppBlaze::blaze_wrap(sv);
   };

   template<typename MT, bool SO>
   SEXP wrap(const blaze::SparseMatrix<MT, SO>& sm) {
   return RcppBlaze::blaze_wrap(sm);
   };
   */
  template<typename Type, size_t N, bool TF>
  SEXP wrap(const blaze::StaticVector<Type, N, TF>& sv) {
    return RcppBlaze::blaze_wrap(sv);
  };

  template<typename Type, size_t N, bool TF>
  SEXP wrap(const blaze::HybridVector<Type, N, TF>& hv) {
    return RcppBlaze::blaze_wrap(hv);
  };

  template<typename Type, bool TF>
  SEXP wrap(const blaze::DynamicVector<Type, TF>& dv) {
    return RcppBlaze::blaze_wrap(dv);
  };

  template<typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool TF>
  SEXP wrap(const blaze::CustomVector<Type, AF, PF, TF>& cv) {
    return RcppBlaze::blaze_wrap(cv);
  };
  /*
   template<typename Type, bool TF>
   SEXP wrap(const blaze::CompressedVector<Type, TF>& cv) {
   return RcppBlaze::blaze_wrap(cv);
   };
   */
  template<typename Type, size_t M, size_t N, bool SO>
  SEXP wrap(const blaze::StaticMatrix<Type, M, N, SO>& sm) {
    return RcppBlaze::blaze_wrap(sm);
  };

  template<typename Type, size_t M, size_t N, bool SO>
  SEXP wrap(const blaze::HybridMatrix<Type, M, N, SO>& hm) {
    return RcppBlaze::blaze_wrap(hm);
  };

  template<typename Type, bool SO>
  SEXP wrap(const blaze::DynamicMatrix<Type, SO>& dm) {
    return RcppBlaze::blaze_wrap(dm);
  };

  template<typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool SO>
  SEXP wrap(const blaze::CustomMatrix<Type, AF, PF, SO>& cm) {
    return RcppBlaze::blaze_wrap(cm);
  };
  /*
   template<typename Type, bool SO>
   SEXP wrap(const blaze::CompressedMatrix<Type, SO>& cm) {
   return RcppBlaze::blaze_wrap(cm);
   };

   template<typename MT, bool SO>
   SEXP wrap(const blaze::IdentityMatrix<MT, SO>& sim) {
   return RcppBlaze::blaze_wrap(sim);
   };

   template<typename MT, bool SO>
   SEXP wrap(const blaze::ZeroMatrix<MT, SO>& szm) {
   return RcppBlaze::blaze_wrap(szm);
   };

   template<typename MT, bool SO>
   SEXP wrap(const blaze::UniformMatrix<MT, SO>& sum) {
   return RcppBlaze::blaze_wrap(sum);
   };
   */

  /*
  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::DiagonalMatrix<MT, SO, DF>& dm) {
    return RcppBlaze::blaze_wrap(dm);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::LowerMatrix<MT, SO, DF>& lm) {
    return RcppBlaze::blaze_wrap(lm);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::UpperMatrix<MT, SO, DF>& um) {
    return RcppBlaze::blaze_wrap(um);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::HermitianMatrix<MT, SO, DF>& hm) {
    return RcppBlaze::blaze_wrap(hm);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::StrictlyLowerMatrix<MT, SO, DF>& slm) {
    return RcppBlaze::blaze_wrap(slm);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::StrictlyUpperMatrix<MT, SO, DF>& sum) {
    return RcppBlaze::blaze_wrap(sum);
  };

  template<typename MT, bool SO, bool DF, bool NF>
  SEXP wrap(const blaze::SymmetricMatrix<MT, SO, DF, NF>& sm) {
    return RcppBlaze::blaze_wrap(sm);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::UniLowerMatrix<MT, SO, DF>& ulm) {
    return RcppBlaze::blaze_wrap(ulm);
  };

  template<typename MT, bool SO, bool DF>
  SEXP wrap(const blaze::UniUpperMatrix<MT, SO, DF>& uum) {
    return RcppBlaze::blaze_wrap(uum);
  };

  // Operations of Vector to Vector
  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::DVecDVecAddExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::DVecDVecCrossExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::DVecDVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::DVecDVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::DVecSVecAddExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::DVecSVecCrossExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::DVecSVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::SVecDVecCrossExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::SVecDVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::SVecSVecCrossExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::DVecSVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::SVecDVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::SVecSVecAddExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::SVecSVecMultExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2, bool TF>
  SEXP wrap(const blaze::SVecSVecSubExpr<VT1, VT2, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::DVecTDVecMultExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::DVecTSVecMultExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::SVecTDVecMultExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT1, typename VT2>
  SEXP wrap(const blaze::SVecTSVecMultExpr<VT1, VT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  // Operations of Matrix to Matrix
  template<typename MT1, typename MT2, bool SO>
  SEXP wrap(const blaze::DMatDMatAddExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2, bool SO>
  SEXP wrap(const blaze::DMatDMatSubExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2, bool SO>
  SEXP wrap(const blaze::DMatSMatAddExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2, bool SO>
  SEXP wrap(const blaze::DMatSMatSubExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatTDMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatTDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatTDMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatTSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatTSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::DMatTSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2, bool SO>
  SEXP wrap(const blaze::SMatDMatSubExpr<MT1, MT2, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatTDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatTDMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TDMatDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TDMatSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TDMatSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TDMatSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TDMatTDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TDMatTSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatDMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatTDMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatTSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatTSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::SMatTSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatSMatSubExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatTSMatAddExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT1, typename MT2>
  SEXP wrap(const blaze::TSMatTSMatMultExpr<MT1, MT2>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  // Operations of Matrix to Vector
  template<typename MT, typename VT>
  SEXP wrap(const blaze::DMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::DMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::TDMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::TDMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::SMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::SMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::TSMatDVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename VT>
  SEXP wrap(const blaze::TSMatSVecMultExpr<MT, VT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  // Operations of Vector to Matrix
  template<typename VT, typename MT>
  SEXP wrap(const blaze::TDVecDMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TDVecTDMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TSVecDMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TSVecTDMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TDVecSMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TDVecTSMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TSVecSMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename MT>
  SEXP wrap(const blaze::TSVecTSMatMultExpr<VT, MT>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  // Operations of Vectors or Matrices
  template<typename VT, bool TF>
  SEXP wrap(const blaze::DVecAbsExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::DVecConjExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::DVecEvalExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::DVecImagExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::DVecRealExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename ST, bool TF>
  SEXP wrap(const blaze::DVecScalarDivExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename ST, bool TF>
  SEXP wrap(const blaze::DVecScalarMultExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::DVecTransExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::SVecAbsExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::SVecConjExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::SVecEvalExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::SVecImagExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::SVecRealExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename ST, bool TF>
  SEXP wrap(const blaze::SVecScalarDivExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, typename ST, bool TF>
  SEXP wrap(const blaze::SVecScalarMultExpr<VT, ST, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename VT, bool TF>
  SEXP wrap(const blaze::SVecTransExpr<VT, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatAbsExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatConjExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatEvalExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatImagExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatInvExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatRealExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename ST, bool SO>
  SEXP wrap(const blaze::DMatScalarDivExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename ST, bool SO>
  SEXP wrap(const blaze::DMatScalarMultExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatTransExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::DMatTransposer<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::SMatAbsExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::SMatConjExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::SMatEvalExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::SMatImagExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::SMatRealExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename ST, bool SO>
  SEXP wrap(const blaze::SMatScalarDivExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, typename ST, bool SO>
  SEXP wrap(const blaze::SMatScalarMultExpr<MT, ST, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename MT, bool SO>
  SEXP wrap(const blaze::SMatTransExpr<MT, SO>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  */
}

#endif
