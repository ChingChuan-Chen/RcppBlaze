// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
// Copyright (C)  2017         Chingchuan Chen
//
// This file is based on RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.
//
// Copyright (C)  2017 Chingchuan Chen
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

#ifndef RcppBlaze__RcppBlazeForward__h
#define RcppBlaze__RcppBlazeForward__h

#include <RcppBlazeConfig.h>
#define STRICT_R_HEADERS
#define R_NO_REMAP
#include <Rinternals.h>
#include <R_ext/Boolean.h>
#include <Rcpp/XPtr.h>
#include <RcppCommon.h>
#include <Rconfig.h>
#include <blaze/Blaze.h>
#include <blaze/system/Version.h>

#ifdef RCPPBLAZE_USE_RCPPLAPACKE
#undef BLAZE_BLAS_MODE
#define BLAZE_BLAS_MODE 1
#undef BLAZE_BLAS_IS_PARALLEL
#define BLAZE_BLAS_IS_PARALLEL 1
#endif

/* forward declarations */
namespace Rcpp {
  /* support for wrap */

  template< typename MT, bool SO > SEXP wrap( const blaze::DenseMatrix<MT,SO>& dm );
  template< typename VT, bool TF > SEXP wrap( const blaze::DenseVector<VT,TF>& dv );
  template< typename VT, bool TF > SEXP wrap( const blaze::SparseVector<VT,TF>& sv );
  template< typename MT, bool SO > SEXP wrap( const blaze::SparseMatrix<MT,SO>& sm );

  template< typename Type, size_t N, bool TF > SEXP wrap( const blaze::StaticVector<Type,N,TF>& );
  template< typename Type, size_t N, bool TF > SEXP wrap( const blaze::HybridVector<Type,N,TF>& );
  template< typename Type, bool TF > SEXP wrap( const blaze::DynamicVector<Type,TF>& );
  template< typename Type, bool AF, bool PF, bool TF > SEXP wrap( const blaze::CustomVector<Type,AF,PF,TF>& );
  template< typename Type, bool TF > SEXP wrap( const blaze::CompressedVector<Type,TF>& );

  template< typename Type, size_t M, size_t N, bool SO > SEXP wrap( const blaze::StaticMatrix<Type,M,N,SO>& );
  template< typename Type, size_t M, size_t N, bool SO > SEXP wrap( const blaze::HybridMatrix<Type,M,N,SO>& );
  template< typename Type, bool SO > SEXP wrap( const blaze::DynamicMatrix<Type,SO>& );
  template< typename Type, bool SO > SEXP wrap( const blaze::DynamicMatrix<Type,SO>& );
  template< typename Type, bool AF, bool PF, bool SO > SEXP wrap( const blaze::CustomMatrix<Type,AF,PF,SO>& );
  template< typename Type, bool SO > SEXP wrap( const blaze::CompressedMatrix<Type,SO>& );

  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::DiagonalMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::LowerMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::UpperMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::HermitianMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::StrictlyLowerMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::StrictlyUpperMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF, bool NF > SEXP wrap( const blaze::SymmetricMatrix<MT,SO,DF,NF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::UniLowerMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::UniUpperMatrix<MT,SO,DF>& );

  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::DenseColumn<MT,SO,SF>& );
  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::DenseRow<MT,SO,SF>& );
  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::DenseSubmatrix<MT,SO,SF>& );
  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::DenseSubvector<MT,SO,SF>& );

  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::SparseColumn<MT,SO,SF>& );
  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::SparseRow<MT,SO,SF>& );
  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::SparseSubmatrix<MT,SO,SF>& );
  template< typename MT, bool SO, bool SF > SEXP wrap( const blaze::SparseSubvector<MT,SO,SF>& );

  namespace traits {

  /* support for as */
  template< typename Type, bool TF > class Exporter< blaze::DynamicVector<Type,TF> >;
  template< typename Type, size_t N, bool TF > class Exporter< blaze::StaticVector<Type,N,TF> >;
  template< typename Type, size_t N, bool TF > class Exporter< blaze::HybridVector<Type,N,TF> >;
  template< typename Type, bool TF > class Exporter< blaze::CustomVector<Type,blaze::unaligned,blaze::unpadded,TF> >;
  template< typename Type, bool TF > class Exporter< blaze::CustomVector<Type,blaze::aligned,blaze::unpadded,TF> >;
  template< typename Type, bool TF > class Exporter< blaze::CustomVector<Type,blaze::unaligned,blaze::padded,TF> >;
  template< typename Type, bool TF > class Exporter< blaze::CustomVector<Type,blaze::aligned,blaze::padded,TF> >;
  template< typename Type > class Exporter< blaze::CompressedVector<Type,blaze::rowVector> >;
  template< typename Type > class Exporter< blaze::CompressedVector<Type,blaze::columnVector> >;

  template< typename Type, bool SO > class Exporter< blaze::DynamicMatrix<Type,SO> >;
  template< typename Type, size_t M, size_t N, bool SO > class Exporter< blaze::StaticMatrix<Type,M,N,SO> >;
  template< typename Type, size_t M, size_t N, bool SO > class Exporter< blaze::HybridMatrix<Type,M,N,SO> >;
  template< typename Type, bool SO > class Exporter< blaze::CustomMatrix<Type,blaze::unaligned,blaze::unpadded,SO> >;
  template< typename Type, bool SO > class Exporter< blaze::CustomMatrix<Type,blaze::aligned,blaze::unpadded,SO> >;
  template< typename Type, bool SO > class Exporter< blaze::CustomMatrix<Type,blaze::unaligned,blaze::padded,SO> >;
  template< typename Type, bool SO > class Exporter< blaze::CustomMatrix<Type,blaze::aligned,blaze::padded,SO> >;
  template< typename Type, bool SO > class Exporter< blaze::CompressedMatrix<Type,SO> >;

  } // namespace traits
}

#endif
