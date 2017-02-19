// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// RcppBlazeWrap.h: Rcpp/Blaze glue
//
// Copyright (C)  2017 Chingchuan Chen
//
// This file is part of RcppBlaze.
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

// most code are modified from the definition of matrix in blaze
namespace Rcpp {

  namespace RcppBlaze {

    /* wrap for the father classes: DenseMatrix, DenseVector, SparseMatrix and SparseVector */

    template< typename MT    // Type of the dense matrix
            , bool SO >      // Transpose flag
    SEXP blaze_wrap( const blaze::DenseMatrix<MT,SO>& dm, ::Rcpp::traits::true_type )
    {
      typedef typename MT::ElementType    ET;
      typedef typename MT::ResultType     RT;
      typedef typename MT::ReturnType     RN;
      typedef typename MT::CompositeType  CT;
      typedef typename blaze::If< blaze::IsExpression<RN>, const RT, CT >::Type  Tmp;
      Tmp A( ~dm );  // Evaluation of the dense matrix operand

      const int RTYPE = ::Rcpp::traits::r_sexptype_traits< ET >::rtype;
      ::Rcpp::Matrix<RTYPE> out( A.rows(), A.columns() );

      if( SO == blaze::rowMajor ) {
        for( size_t i=0UL; i<A.rows(); ++i ) {
          for( size_t j=0UL; j<A.columns(); ++j ) {
            out(i, j) = ::Rcpp::internal::caster< ET, typename Rcpp::traits::storage_type<RTYPE>::type >( A(i, j) );
          }
        }
      } else {
        for( size_t j=0UL; j<A.columns(); ++j ) {
          for( size_t i=0UL; i<A.rows(); ++i ) {
            out(i, j) = ::Rcpp::internal::caster< ET, typename Rcpp::traits::storage_type<RTYPE>::type >( A(i, j) );
          }
        }
      }
      return out;
    }

    template< typename MT    // Type of the dense matrix
            , bool SO >      // Transpose flag
    SEXP blaze_wrap( const blaze::DenseMatrix<MT,SO>& dm, ::Rcpp::traits::false_type )
    {
      typedef typename MT::ElementType    ET;
      typedef typename MT::ResultType     RT;
      typedef typename MT::ReturnType     RN;
      typedef typename MT::CompositeType  CT;
      typedef typename blaze::If< blaze::IsExpression<RN>, const RT, CT >::Type  Tmp;
      Tmp A( ~dm );  // Evaluation of the dense matrix operand

      const int RTYPE = ::Rcpp::traits::r_sexptype_traits< ET >::rtype;
      ::Rcpp::Matrix<RTYPE> out( A.rows(), A.columns() );

      if( SO == blaze::rowMajor ) {
        for( size_t i=0UL; i<A.rows(); ++i ) {
          for( size_t j=0UL; j<A.columns(); ++j ) {
            out(i, j) = A(i, j);
          }
        }
      } else {
        for( size_t j=0UL; j<A.columns(); ++j ) {
          for( size_t i=0UL; i<A.rows(); ++i ) {
            out(i, j) = A(i, j);
          }
        }
      }
      return out;
    }

    template< typename MT    // Type of the dense matrix
      , bool SO >      // Transpose flag
    SEXP blaze_wrap( const blaze::DenseMatrix<MT,SO>& dm )
    {
      typedef typename MT::ElementType ET;
      return blaze_wrap<MT,SO>( dm, typename ::Rcpp::traits::r_sexptype_needscast< ET >() );
    }

    template< typename VT  // Type of the dense vector
            , bool TF >    // Transpose flag
    SEXP blaze_wrap( const blaze::DenseVector<VT,TF>& dv, ::Rcpp::traits::true_type )
    {
      typedef typename VT::CompositeType  CT;
      typedef typename VT::ElementType ET;
      CT a( ~dv );  // Evaluation of the dense vector operand

      const int RTYPE = ::Rcpp::traits::r_sexptype_traits< ET >::rtype;
      ::Rcpp::Vector<RTYPE> out( a.size() );

      for( size_t i=0UL; i<a.size(); ++i )
        out[i] = ::Rcpp::internal::caster< ET, typename Rcpp::traits::storage_type<RTYPE>::type >( a[i] );
      return out;
    }

    template< typename VT  // Type of the dense vector
            , bool TF >    // Transpose flag
    SEXP blaze_wrap( const blaze::DenseVector<VT,TF>& dv, ::Rcpp::traits::false_type )
    {
      typedef typename VT::CompositeType  CT;
      typedef typename VT::ElementType ET;
      CT a( ~dv );  // Evaluation of the dense vector operand

      const int RTYPE = ::Rcpp::traits::r_sexptype_traits< ET >::rtype;
      ::Rcpp::Vector<RTYPE> out( a.size() );

      for( size_t i=0UL; i<a.size(); ++i )
        out[i] = a[i];
      return out;
    }

    template< typename VT  // Type of the dense vector
      , bool TF >    // Transpose flag
    SEXP blaze_wrap( const blaze::DenseVector<VT,TF>& dv )
    {
      typedef typename VT::ElementType ET;
      return blaze_wrap<VT,TF>( dv, typename ::Rcpp::traits::r_sexptype_needscast< ET >() );
    }

    template< typename MT  // Type of the sparse matrix
            , bool SO >    // Storage order
    SEXP blaze_wrap( const blaze::SparseMatrix<MT,SO>& sm )
    {
      typedef typename MT::ResultType     RT;
      typedef typename MT::ReturnType     RN;
      typedef typename MT::CompositeType  CT;
      typedef typename blaze::If< blaze::IsExpression<RN>, const RT, CT >::Type  Tmp;
      typedef typename blaze::RemoveReference<Tmp>::Type::ConstIterator   ConstIterator;

      Tmp A( ~sm );  // Evaluation of the sparse matrix operand

      const int RTYPE = ::Rcpp::traits::r_sexptype_traits<typename MT::ElementType>::rtype;
      std::string klass ;
      switch( RTYPE ){
      case REALSXP: klass = (SO == blaze::rowMajor) ? "dgRMatrix" : "dgCMatrix" ; break ;
      // case INTSXP: klass = (SO == blaze::rowMajor) ? "igRMatrix" : "igCMatrix" ; break ; // class not exported
      case LGLSXP: klass = (SO == blaze::rowMajor) ? "lgRMatrix" : "lgCMatrix" ; break ;
      default:
        throw std::invalid_argument( "RTYPE not matched in conversion to sparse matrix" ) ;
      }

      const size_t index( ( SO == blaze::rowMajor )?( A.rows() ):( A.columns() ) );

      ::Rcpp::Vector<RTYPE> x( A.nonZeros() );
      ::Rcpp::IntegerVector idx( A.nonZeros() );
      ::Rcpp::IntegerVector p( index + 1UL );

      size_t k = 0UL;
      for( size_t i=0UL; i<index; ++i ) {
        p[i] = (int) k;
        for( ConstIterator element=A.begin(i); element!=A.end(i); ++element ) {
          x[k] = element->value();
          idx[k] = element->index();
          ++k;
        }
      }
      p[index] = (int) A.nonZeros();

      ::Rcpp::S4 s(klass);
      s.slot("Dim") = ::Rcpp::Dimension(A.rows(), A.columns());
      s.slot( (SO == blaze::rowMajor) ? "j" : "i" ) = idx;
      s.slot("p") = p;
      s.slot("x") = x;
      return s;
    }

    template< typename VT  // Type of the sparse vector
            , bool TF >    // Transpose flag
    SEXP blaze_wrap( const blaze::SparseVector<VT,TF>& sv )
    {
      typedef typename VT::ElementType    ET;
      typedef typename VT::CompositeType  CT;
      typedef typename blaze::RemoveReference<CT>::Type::ConstIterator  ConstIterator;

      CT a( ~sv );  // Evaluation of the sparse vector operand

      const int RTYPE = ::Rcpp::traits::r_sexptype_traits<typename VT::ElementType>::rtype;
      std::string klass ;
      switch( RTYPE ){
      case REALSXP: klass = "dgCMatrix" ; break ;
      // case INTSXP: klass = (SO == blaze::rowMajor) ? "igRMatrix" : "igCMatrix" ; break ; // class not exported
      case LGLSXP: klass = "lgCMatrix" ; break ;
      default:
        throw std::invalid_argument( "RTYPE not matched in conversion to sparse matrix" ) ;
      }

      ::Rcpp::Vector<RTYPE> x( a.nonZeros() );
      ::Rcpp::IntegerVector idx( a.nonZeros() );
      size_t pSize = 2UL;
      if (TF == blaze::rowVector)
        pSize = a.size() + 1UL;
      ::Rcpp::IntegerVector p( pSize, 0 );

      size_t k = 0UL;
      if (TF == blaze::columnVector) {
        p[1UL] = a.nonZeros();
        for( ConstIterator element=a.begin(); element!=a.end(); ++element ) {
          x[k] = element->value();
          idx[k] = element->index();
          ++k;
        }
      } else if (TF == blaze::rowVector) {
        for( ConstIterator element=a.begin(); element!=a.end(); ++element ) {
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
  } // namespace RcppBlaze

  /* wrap for the child classes of DenseMatrix, DenseVector, SparseMatrix and SparseVector */

  template< typename MT, bool SO >
  SEXP wrap( const blaze::DenseMatrix<MT,SO>& dm ){
    return RcppBlaze::blaze_wrap( dm ) ;
  };

  template< typename VT, bool TF >
  SEXP wrap( const blaze::DenseVector<VT,TF>& dv ){
    return RcppBlaze::blaze_wrap( dv ) ;
  };

  template< typename VT, bool TF >
  SEXP wrap( const blaze::SparseVector<VT,TF>& sv ){
    return RcppBlaze::blaze_wrap( sv ) ;
  };

  template< typename MT, bool SO >
  SEXP wrap( const blaze::SparseMatrix<MT,SO>& sm ){
    return RcppBlaze::blaze_wrap( sm ) ;
  };

  template< typename Type, size_t N, bool TF >
  SEXP wrap( const blaze::StaticVector<Type,N,TF>& sv ){
    return RcppBlaze::blaze_wrap( sv ) ;
  };

  template< typename Type, size_t N, bool TF >
  SEXP wrap( const blaze::HybridVector<Type,N,TF>& hv ){
    return RcppBlaze::blaze_wrap( hv ) ;
  };

  template< typename Type, bool TF >
  SEXP wrap( const blaze::DynamicVector<Type,TF>& dv ){
    return RcppBlaze::blaze_wrap( dv ) ;
  };

  template< typename Type, bool AF, bool PF, bool TF >
  SEXP wrap( const blaze::CustomVector<Type,AF,PF,TF>& cv ){
    return RcppBlaze::blaze_wrap( cv ) ;
  };

  template< typename Type, bool TF >
  SEXP wrap( const blaze::CompressedVector<Type,TF>& cv ){
    return RcppBlaze::blaze_wrap( cv ) ;
  };

  template< typename Type, size_t M, size_t N, bool SO >
  SEXP wrap( const blaze::StaticMatrix<Type,M,N,SO>& sm){
    return RcppBlaze::blaze_wrap( sm ) ;
  };

  template< typename Type, size_t M, size_t N, bool SO >
  SEXP wrap( const blaze::HybridMatrix<Type,M,N,SO>& hm ){
    return RcppBlaze::blaze_wrap( hm ) ;
  };

  template< typename Type, bool SO >
  SEXP wrap( const blaze::DynamicMatrix<Type,SO>& dm ){
    return RcppBlaze::blaze_wrap( dm ) ;
  };

  template< typename Type, bool AF, bool PF, bool SO >
  SEXP wrap( const blaze::CustomMatrix<Type,AF,PF,SO>& cm ){
    return RcppBlaze::blaze_wrap( cm ) ;
  };

  template< typename Type, bool SO >
  SEXP wrap( const blaze::CompressedMatrix<Type,SO>& cm ){
    return RcppBlaze::blaze_wrap( cm ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::DenseColumn<MT,SO,SF>& dc ){
    return RcppBlaze::blaze_wrap( dc ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::DenseRow<MT,SO,SF>& dr ){
    return RcppBlaze::blaze_wrap( dr ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::DenseSubmatrix<MT,SO,SF>& dsm ){
    return RcppBlaze::blaze_wrap( dsm ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::DenseSubvector<MT,SO,SF>& dsv ){
    return RcppBlaze::blaze_wrap( dsv ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::SparseColumn<MT,SO,SF>& sc ){
    return RcppBlaze::blaze_wrap( sc ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::SparseRow<MT,SO,SF>& sr ){
    return RcppBlaze::blaze_wrap( sr ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::SparseSubmatrix<MT,SO,SF>& ssm ){
    return RcppBlaze::blaze_wrap( ssm ) ;
  };

  template< typename MT, bool SO, bool SF >
  SEXP wrap( const blaze::SparseSubvector<MT,SO,SF>& ssv ){
    return RcppBlaze::blaze_wrap( ssv ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::DiagonalMatrix<MT,SO,DF>& dm ){
    return RcppBlaze::blaze_wrap( dm ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::LowerMatrix<MT,SO,DF>& lm ){
    return RcppBlaze::blaze_wrap( lm ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::UpperMatrix<MT,SO,DF>& um ){
    return RcppBlaze::blaze_wrap( um ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::HermitianMatrix<MT,SO,DF>& hm ){
    return RcppBlaze::blaze_wrap( hm ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::StrictlyLowerMatrix<MT,SO,DF>& slm ){
    return RcppBlaze::blaze_wrap( slm ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::StrictlyUpperMatrix<MT,SO,DF>& sum ){
    return RcppBlaze::blaze_wrap( sum ) ;
  };

  template< typename MT, bool SO, bool DF, bool NF >
  SEXP wrap( const blaze::SymmetricMatrix<MT,SO,DF,NF>& sm ){
    return RcppBlaze::blaze_wrap( sm ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::UniLowerMatrix<MT,SO,DF>& ulm ){
    return RcppBlaze::blaze_wrap( ulm ) ;
  };

  template< typename MT, bool SO, bool DF >
  SEXP wrap( const blaze::UniUpperMatrix<MT,SO,DF>& uum ){
    return RcppBlaze::blaze_wrap( uum ) ;
  };
}

#endif
