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

#include "RcppBlazeForward.h"
#include <Rcpp.h>

// most code are modified from the definition of matrix in blaze

// TODO:
// Convert IdentityMatrix, ZeroMatrix
// Convert Symmetric Matrices, Hermitian Matrices, Triangular Matrices
// Convert Submatrix, Subvector, Row, Column, Rows, Columns, Band

namespace Rcpp {

  namespace RcppBlaze {

    template<typename VT, bool TF>
    SEXP blaze_wrap(const blaze::DenseVector<VT, TF>& x) {
      typedef typename VT::ElementType ET;
      const int RTYPE = Rcpp::traits::r_sexptype_traits<ET>::rtype;
      const ET* data_pointer = blaze::data(*x);
      Rcpp::Vector<RTYPE> out = Rcpp::wrap(data_pointer, data_pointer + (*x).size());
      return out;
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



  template <typename Type, bool TF>
  SEXP wrap(const blaze::DynamicVector<Type, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename Type, size_t N, bool TF>
  SEXP wrap(const blaze::StaticVector<Type, N, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename Type, size_t N, bool TF>
  SEXP wrap(const blaze::HybridVector<Type, N, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  template<typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool TF>
  SEXP wrap(const blaze::CustomVector<Type, AF, PF, TF>& x) {
    return RcppBlaze::blaze_wrap(x);
  };

  /*
   template<typename Type, bool TF>
   SEXP wrap(const blaze::CompressedVector<Type, TF>& cv) {
   return RcppBlaze::blaze_wrap(cv);
   };

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
   */
}

#endif
