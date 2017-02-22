// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
// Copyright (C)  2017 Chingchuan Chen
//
// This file is based on RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.
//
// RcppBlazeAs.h: Rcpp/Blaze glue
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

#ifndef RcppBlaze__RcppBlazeAs__h
#define RcppBlaze__RcppBlazeAs__h

// most code are modified from the Exporter class definition in Rcpp / RcppArmadillo / RcppEigen
namespace Rcpp {

  namespace traits {

    // Provides only blaze::DynamicVector<Type,TF> export
    template< typename Type, bool TF >
    class Exporter< blaze::DynamicVector<Type,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {}
      blaze::DynamicVector<Type,TF> get() {
        blaze::DynamicVector<Type,TF> result( (size_t) vec.size() );
        for( size_t i=0UL; i<(size_t)vec.size(); ++i )
          result[i] = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( vec[i] );
        return result ;
      }
    };

    // Provides only blaze::StaticVector<Type,N,TF> export
    template< typename Type, size_t N, bool TF >
    class Exporter< blaze::StaticVector<Type,N,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {
        if( vec.size() != N )
          throw std::invalid_argument("Size is not match.");
      }
      blaze::StaticVector<Type,N,TF> get() {
        blaze::StaticVector<Type,N,TF> result;
        for( size_t i=0UL; i<(size_t)vec.size(); ++i )
          result[i] = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( vec[i] );
        return result ;
      }
    };

    // Provides only blaze::HybridVector<Type,N,TF> export
    template< typename Type, size_t N, bool TF >
    class Exporter< blaze::HybridVector<Type,N,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {
        if( vec.size() != N )
          throw std::invalid_argument("Size is not match.");
      }
      blaze::HybridVector<Type,N,TF> get() {
        blaze::HybridVector<Type,N,TF> result( (size_t) vec.size() );
        for( size_t i=0UL; i<(size_t)vec.size(); ++i )
          result[i] = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( vec[i] );
        return result ;
      }
    };

    // Provides only blaze::CustomVector<Type,blaze::unaligned,blaze::unpadded,TF> export
    template< typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type,blaze::unaligned,blaze::unpadded,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type,blaze::unaligned,blaze::unpadded,TF> get() {
        blaze::CustomVector<Type,blaze::unaligned,blaze::unpadded,TF> result( reinterpret_cast<Type*>(vec.begin()), (size_t) vec.size() );
        return result;
      }
    };

    // Provides only blaze::CustomVector<Type,blaze::aligned, blaze::unpadded,TF> export
    template< typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type,blaze::aligned, blaze::unpadded,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type,blaze::aligned,blaze::unpadded,TF> get() {
        blaze::CustomVector<Type,blaze::aligned,blaze::unpadded,TF> result(
            blaze::allocate<Type>( (size_t)vec.size() ), (size_t) vec.size(), blaze::Deallocate() );
        for( size_t i=0UL; i<(size_t)vec.size(); ++i )
          result[i] = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( vec[i] );
        return result;
      }
    };

    // Provides only blaze::CustomVector<Type,blaze::unaligned, blaze::padded,TF> export
    template< typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type,blaze::unaligned, blaze::padded,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type,blaze::unaligned,blaze::padded,TF> get() {
        size_t paddedSize = (size_t) vec.size() + blaze::IntrinsicTrait<Type>::size -
          (size_t) vec.size() % blaze::IntrinsicTrait<Type>::size;
        blaze::CustomVector<Type,blaze::unaligned,blaze::padded,TF> result(
            new Type[paddedSize], (size_t) vec.size(), paddedSize, blaze::ArrayDelete() );
        for( size_t i=0UL; i<(size_t)vec.size(); ++i )
          result[i] = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( vec[i] );
        return result;
      }
    };

    // Provides only blaze::CustomVector<Type,blaze::aligned, blaze::padded,TF> export
    template< typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type,blaze::aligned, blaze::padded,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Vector<RTYPE> vec ;

    public:
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type,blaze::aligned,blaze::padded,TF> get() {
        size_t paddedSize = (size_t) vec.size() + blaze::IntrinsicTrait<Type>::size -
          (size_t) vec.size() % blaze::IntrinsicTrait<Type>::size;
        blaze::CustomVector<Type,blaze::aligned,blaze::padded,TF> result(
            blaze::allocate<Type>(paddedSize), (size_t) vec.size(), paddedSize, blaze::Deallocate() );
        for( size_t i=0UL; i<(size_t)vec.size(); ++i )
          result[i] = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( vec[i] );
        return result;
      }
    };

    // Provides only blaze::CompressedVector<Type,blaze::rowVector> export
    template< typename Type >
    class Exporter< blaze::CompressedVector<Type,blaze::rowVector> > {
      ::Rcpp::S4 vec;

    public:
      Exporter(SEXP x) : vec(x) {
        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        if( dims[0] != 1 )
          throw std::invalid_argument("Not a sparse row vector.");
      }
      blaze::CompressedVector<Type,blaze::rowVector> get() {
        const int  RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        ::Rcpp::Vector<RTYPE> x = vec.slot("x");

        blaze::CompressedVector<Type,blaze::rowVector> result( (size_t) dims[1], (size_t) x.size() );
        ::Rcpp::IntegerVector p = vec.slot("p");
        size_t i = 0UL;
        for( size_t k=1UL; k<(size_t)p.size()-1UL; ++k ) {
          if (p[k+1] - p[k] > 0) {
            result[ (size_t) k] = x[i];
            ++i;
          }
        }
        return result;
      }
    };

    // Provides only blaze::CompressedVector<Type,blaze::columnVector> export
    template< typename Type >
    class Exporter< blaze::CompressedVector<Type,blaze::columnVector> > {
      ::Rcpp::S4 vec;

    public:
      Exporter(SEXP x) : vec(x) {
        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        if( dims[1] != 1 )
          throw std::invalid_argument("Not a sparse column vector.");
      }
      blaze::CompressedVector<Type,blaze::columnVector> get() {
        const int  RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        ::Rcpp::Vector<RTYPE> x = vec.slot("x");

        blaze::CompressedVector<Type,blaze::columnVector> result( (size_t) dims[0], (size_t) x.size() );
        ::Rcpp::IntegerVector i = vec.slot("i");
        for( size_t k=0UL; k<(size_t)x.size(); ++k )
          result[ (size_t) i[k] ] = x[k];
        return result;
      }
    };

    // Provides only blaze::DynamicMatrix<Type,SO> export
    template< typename Type, bool SO >
    class Exporter< blaze::DynamicMatrix<Type,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::DynamicMatrix<Type,SO> get() {
        blaze::DynamicMatrix<Type,SO> result( (size_t) mat.nrow(), (size_t) mat.ncol() );
        if( SO == blaze::rowMajor ) {
          for( size_t i=0UL; i<result.rows(); ++i ) {
            for( size_t j=0UL; j<result.columns(); ++j )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
        } else {
          for( size_t j=0UL; j<result.columns(); ++j ) {
            for( size_t i=0UL; i<result.rows(); ++i )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
        }
        return result ;
      }
    };

    // Provides only blaze::StaticMatrix<Type,M,N,SO> export
    template< typename Type, size_t M, size_t N, bool SO >
    class Exporter< blaze::StaticMatrix<Type,M,N,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {
        if ( (size_t) mat.nrow() != M && (size_t) mat.ncol() != N )
          throw std::invalid_argument("Dimension of matrix is not match.");
      }
      blaze::StaticMatrix<Type,M,N,SO> get() {
        blaze::StaticMatrix<Type,M,N,SO> result;
        if( SO == blaze::rowMajor ) {
          for( size_t i=0UL; i<result.rows(); ++i ) {
            for( size_t j=0UL; j<result.columns(); ++j )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
        } else {
          for( size_t j=0UL; j<result.columns(); ++j ) {
            for( size_t i=0UL; i<result.rows(); ++i )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
        }
        return result ;
      }
    };

    // Provides only blaze::HybridMatrix<Type,M,N,SO> export
    template< typename Type, size_t M, size_t N, bool SO >
    class Exporter< blaze::HybridMatrix<Type,M,N,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {
        if ( (size_t) mat.nrow() != M && (size_t) mat.ncol() != N )
          throw std::invalid_argument("Dimension of matrix is not match.");
      }
      blaze::HybridMatrix<Type,M,N,SO> get() {
        blaze::HybridMatrix<Type,M,N,SO> result( (size_t) mat.nrow(), (size_t) mat.ncol() );
        if( SO == blaze::rowMajor ) {
          for( size_t i=0UL; i<result.rows(); ++i ) {
            for( size_t j=0UL; j<result.columns(); ++j )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
        } else {
          for( size_t j=0UL; j<result.columns(); ++j ) {
            for( size_t i=0UL; i<result.rows(); ++i )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
        }
        return result ;
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> export
    template< typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> get() {
        blaze::CustomMatrix<Type,blaze::unaligned,blaze::unpadded,SO> result( reinterpret_cast<Type*>(mat.begin()), (size_t) mat.nrow(), (size_t) mat.ncol());
        return result ;
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> export
    template< typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> get() {
        if( SO == blaze::rowMajor ) {
          size_t actualCols = (size_t) mat.ncol() + blaze::IntrinsicTrait<Type>::size -
            (size_t) mat.ncol() % blaze::IntrinsicTrait<Type>::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> result(
              blaze::allocate<Type>( actualCols * (size_t)mat.nrow() ), (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualCols, blaze::Deallocate() );
          for( size_t i=0UL; i<result.rows(); ++i ) {
            for( size_t j=0UL; j<result.columns(); ++j )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
          return result;
        } else {
          size_t actualRows = (size_t) mat.nrow() + blaze::IntrinsicTrait<Type>::size -
            (size_t) mat.nrow() % blaze::IntrinsicTrait<Type>::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> result(
              blaze::allocate<Type>( actualRows * (size_t)mat.ncol() ), (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualRows, blaze::Deallocate() );
          for( size_t j=0UL; j<result.columns(); ++j ) {
            for( size_t i=0UL; i<result.rows(); ++i )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
          return result;
        }
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> export
    template< typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> get() {
        if( SO == blaze::rowMajor ) {
          size_t actualCols = (size_t) mat.ncol() + blaze::IntrinsicTrait<Type>::size -
            (size_t) mat.ncol() % blaze::IntrinsicTrait<Type>::size;
          blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> result(
              new Type[ actualCols * (size_t)mat.nrow() ], (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualCols, blaze::ArrayDelete() );
          for( size_t i=0UL; i<result.rows(); ++i ) {
            for( size_t j=0UL; j<result.columns(); ++j )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
          return result;
        } else {
          size_t actualRows = (size_t) mat.nrow() + blaze::IntrinsicTrait<Type>::size -
            (size_t) mat.nrow() % blaze::IntrinsicTrait<Type>::size;
          blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> result(
              new Type[ actualRows * (size_t)mat.ncol() ], (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualRows, blaze::ArrayDelete() );
          for( size_t j=0UL; j<result.columns(); ++j ) {
            for( size_t i=0UL; i<result.rows(); ++i )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
          return result;
        }
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> export
    template< typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype ;
      ::Rcpp::Matrix<RTYPE> mat ;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> get() {
        if( SO == blaze::rowMajor ) {
          size_t actualCols = (size_t) mat.ncol() + blaze::IntrinsicTrait<Type>::size -
            (size_t) mat.ncol() % blaze::IntrinsicTrait<Type>::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> result(
              blaze::allocate<Type>( actualCols * (size_t)mat.nrow() ), (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualCols, blaze::Deallocate() );
          for( size_t i=0UL; i<result.rows(); ++i ) {
            for( size_t j=0UL; j<result.columns(); ++j )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
          return result;
        } else {
          size_t actualRows = (size_t) mat.nrow() + blaze::IntrinsicTrait<Type>::size -
            (size_t) mat.nrow() % blaze::IntrinsicTrait<Type>::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> result(
              blaze::allocate<Type>( actualRows * (size_t)mat.ncol() ), (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualRows, blaze::Deallocate() );
          for( size_t j=0UL; j<result.columns(); ++j ) {
            for( size_t i=0UL; i<result.rows(); ++i )
              result(i ,j) = ::Rcpp::internal::caster< typename Rcpp::traits::storage_type<RTYPE>::type, Type >( mat(i ,j) );
          }
          return result;
        }
      }
    };

    // Provides only blaze::CompressedMatrix<Type,SO> export
    template< typename Type, bool SO >
    class Exporter< blaze::CompressedMatrix<Type,SO> > {
      ::Rcpp::S4 mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CompressedMatrix<Type,SO> get() {
        const int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
        ::Rcpp::IntegerVector dims = mat.slot("Dim");
        ::Rcpp::Vector<RTYPE> x = mat.slot("x");
        ::Rcpp::IntegerVector i = mat.slot("i");
        ::Rcpp::IntegerVector p = mat.slot("p");
        blaze::CompressedMatrix<Type,SO> result( (size_t) dims[0], (size_t) dims[1] );

        int j = 0;
        for( size_t k=0UL; k<(size_t)x.size(); ++k ) {
          if ( k == (size_t) p[j+1] )
            ++j;
          result( i[k], j ) = x[k];
        }
        return result;
      }
    };
  } // end traits
}

#endif
