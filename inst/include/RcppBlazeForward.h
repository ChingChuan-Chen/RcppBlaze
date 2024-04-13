// Copyright (C)  2010 - 2024 Dirk Eddelbuettel and Romain Francois
// Copyright (C)  2017 - 2024 Ching-Chuan Chen
//
// This file is based on files from RcppArmadillo.
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.

#ifndef RcppBlaze__RcppBlazeForward__h
#define RcppBlaze__RcppBlazeForward__h

#include <RcppBlazeConfig.h>

#ifndef STRICT_R_HEADERS
#define STRICT_R_HEADERS
#endif

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#include <Rinternals.h>
#include <R_ext/Boolean.h>
#include <Rcpp/XPtr.h>
#include <Rconfig.h>

#include <blaze/Blaze.h>

/* forward declarations */
namespace Rcpp {

  /* support wrapping blaze dense vectors */
  template <typename Type, bool TF> SEXP wrap(const blaze::DynamicVector<Type, TF>&);
  template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::StaticVector<Type, N, TF, AF, PF>&);
  template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::HybridVector<Type, N, TF, AF, PF>&);
  template <typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool TF>
  SEXP wrap(const blaze::CustomVector<Type, AF, PF, TF>&);
  template <typename Type, bool TF> SEXP wrap(const blaze::UniformVector<Type, TF>&);

  /* support wrapping blaze dense matrices */
  template <typename Type, bool SO> SEXP wrap(const blaze::DynamicMatrix<Type, SO>&);
  template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::StaticMatrix<Type, M, N, SO, AF, PF>&);
  template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
  SEXP wrap(const blaze::HybridMatrix<Type, M, N, SO, AF, PF>&);
  template <typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool SO>
  SEXP wrap(const blaze::CustomMatrix<Type, AF, PF, SO>&);
  template <typename Type, bool SO> SEXP wrap(const blaze::UniformMatrix<Type, SO>&);

  /* support wrapping blaze sparse vectors */
  template <typename Type, bool TF> SEXP wrap(const blaze::CompressedVector<Type,TF>&);
  template <typename Type, bool TF> SEXP wrap(const blaze::ZeroVector<Type,TF>&);

  /* support wrapping blaze sparse matrices */
  template <typename Type, bool SO> SEXP wrap(const blaze::CompressedMatrix<Type,SO>&);
  template <typename Type, bool SO> SEXP wrap(const blaze::IdentityMatrix<Type,SO>&);
  template <typename Type, bool SO> SEXP wrap(const blaze::ZeroMatrix<Type,SO>&);
  /*
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::DiagonalMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::LowerMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::UpperMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::HermitianMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::StrictlyLowerMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::StrictlyUpperMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF, bool NF > SEXP wrap( const blaze::SymmetricMatrix<MT,SO,DF,NF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::UniLowerMatrix<MT,SO,DF>& );
  template< typename MT, bool SO, bool DF > SEXP wrap( const blaze::UniUpperMatrix<MT,SO,DF>& );
  */

  namespace traits {

    /* support for as for blaze dense vectors */
    template <typename Type, bool TF> class Exporter< blaze::DynamicVector<Type, TF> >;
    template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::HybridVector<Type, N, TF, AF, PF> >;
    template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticVector<Type, N, TF, AF, PF> >;
    template <typename Type, bool TF> class Exporter< blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> >;
    template <typename Type, bool TF> class Exporter< blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> >;
    template <typename Type, bool TF> class Exporter< blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> >;
    template <typename Type, bool TF> class Exporter< blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> >;

    /* support for as for blaze dense matrices */
    template <typename Type, bool SO> class Exporter< blaze::DynamicMatrix<Type, SO> >;
    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::HybridMatrix<Type, M, N, SO, AF, PF> >;
    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticMatrix<Type, M, N, SO, AF, PF> >;
    template <typename Type, bool SO> class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> >;
    template <typename Type, bool SO> class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::unpadded, SO> >;
    template <typename Type, bool SO> class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO> >;
    template <typename Type, bool SO> class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO> >;

    /* support for as for blaze sparse vectors */
    template <typename Type, bool TF> class Exporter< blaze::CompressedVector<Type, TF> >;
    template <typename Type, bool TF> class Exporter< blaze::ZeroVector<Type, TF> >;

    /* support for as for blaze sparse matrices */
    template <typename Type, bool SO> class Exporter< blaze::CompressedMatrix<Type, SO> >;
    template <typename Type, bool SO> class Exporter< blaze::ZeroMatrix<Type, SO> >;

  } // namespace traits
}

#endif
