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

#include <RcppCommon.h>
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

  /* support blaze matrices adaptors */
  template <typename MT, bool SO, bool DF, bool NF > SEXP wrap(const blaze::SymmetricMatrix<MT, SO, DF, NF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::HermitianMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::LowerMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::UniLowerMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::StrictlyLowerMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::UpperMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::UniUpperMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::StrictlyUpperMatrix<MT, SO, DF>&);
  template <typename MT, bool SO, bool DF> SEXP wrap(const blaze::DiagonalMatrix<MT, SO, DF>&);

  /* support blaze views for vectors */
  template <typename VT, blaze::AlignmentFlag AF, bool TF, bool DF, size_t... CSAs>
  SEXP wrap(const blaze::Subvector<VT, AF, TF, DF, CSAs...>&);
  template <typename VT, bool TF, bool DF, typename... CEAs>
  SEXP wrap(const blaze::Elements<VT, TF, DF, CEAs...>&);

  /* support blaze views for matrices */
  template <typename MT, bool SO, bool DF, bool SF, size_t... CRAs>
  SEXP wrap(const blaze::Row<MT, SO, DF, SF, CRAs...>&);
  template <typename MT, bool SO, bool DF, bool SF, size_t... CCAs>
  SEXP wrap(const blaze::Column<MT, SO, DF, SF, CCAs...>&);
  template <typename MT, bool TF, bool DF, bool MF, ptrdiff_t... CBAs>
  SEXP wrap(const blaze::Band<MT, TF, DF, MF, CBAs...>&);
  template <typename MT, blaze::AlignmentFlag AF, bool SO, bool DF, size_t... CSAs>
  SEXP wrap(const blaze::Submatrix<MT, AF, SO, DF, CSAs...>&);
  template <typename MT, bool SO, bool DF, bool SF, typename... CRAs>
  SEXP wrap(const blaze::Rows<MT, SO, DF, SF, CRAs...>&);
  template <typename MT, bool SO, bool DF, bool SF, typename... CCAs>
  SEXP wrap(const blaze::Columns<MT, SO, DF, SF, CCAs...>&);

  /* support blaze Expressions */
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatAddExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatKronExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, typename OP, bool SO> SEXP wrap(const blaze::DMatDMatMapExpr<MT1, MT2, OP, SO>&);
  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::DMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSchurExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSolveExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSubExpr<MT1, MT2, SO>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::DMatDVecMultExpr<MT, VT>&);
  template <typename MT, typename VT, bool TF> SEXP wrap(const blaze::DMatDVecSolveExpr<MT, VT, TF>&);
  template <typename MT, typename OP, bool SO> SEXP wrap(const blaze::DMatMapExpr<MT, OP, SO>&);
  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::DMatScalarDivExpr<MT, ST, SO>&);
  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::DMatScalarMultExpr<MT, ST, SO>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatAddExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatKronExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::DMatSMatMultExpr<MT1, MT2, SF, HF, LF, UF>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::DMatSMatSchurExpr<MT1, MT2>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatSubExpr<MT1, MT2, SO>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::DMatSVecMultExpr<MT, VT>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecAddExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecCrossExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecDivExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecKronExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, typename OP, bool TF> SEXP wrap(const blaze::DVecDVecMapExpr<VT1, VT2, OP, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecMultExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, typename OP> SEXP wrap(const blaze::DVecDVecOuterExpr<VT1, VT2, OP>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecSubExpr<VT1, VT2, TF>&);
  template <typename VT, typename OP, bool TF> SEXP wrap(const blaze::DVecMapExpr<VT, OP, TF>&);
  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::DVecScalarDivExpr<VT, ST, TF>&);
  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::DVecScalarMultExpr<VT, ST, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecAddExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecCrossExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecKronExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecMultExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2> SEXP wrap(const blaze::DVecSVecOuterExpr<VT1, VT2>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecSubExpr<VT1, VT2, TF>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::SMatDMatKronExpr<MT1, MT2, SO>&);
  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::SMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatDMatSchurExpr<MT1, MT2>&);
  template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::SMatDMatSubExpr<MT1, MT2, SO>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::SMatDVecMultExpr<MT, VT>&);
  template <typename MT, typename OP, bool SO> SEXP wrap(const blaze::SMatMapExpr<MT, OP, SO>&);
  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::SMatScalarDivExpr<MT, ST, SO>&);
  template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::SMatScalarMultExpr<MT, ST, SO>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatAddExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatKronExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatMultExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatSchurExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatSubExpr<MT1, MT2>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::SMatSVecMultExpr<MT, VT>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecCrossExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecDivExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecKronExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecMultExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2> SEXP wrap(const blaze::SVecDVecOuterExpr<VT1, VT2>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecSubExpr<VT1, VT2, TF>&);
  template <typename VT, typename OP, bool TF> SEXP wrap(const blaze::SVecMapExpr<VT, OP, TF>&);
  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::SVecScalarDivExpr<VT, ST, TF>&);
  template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::SVecScalarMultExpr<VT, ST, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecAddExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecCrossExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecKronExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecMultExpr<VT1, VT2, TF>&);
  template <typename VT1, typename VT2> SEXP wrap(const blaze::SVecSVecOuterExpr<VT1, VT2>&);
  template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecSubExpr<VT1, VT2, TF>&);
  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TDMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::TDMatDVecMultExpr<MT, VT>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TDMatSMatAddExpr<MT1, MT2>&);
  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TDMatSMatMultExpr<MT1, MT2, SF, HF, LF, UF>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TDMatSMatSubExpr<MT1, MT2>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::TDMatSVecMultExpr<MT, VT>&);
  template <typename VT, typename MT> SEXP wrap(const blaze::TDVecDMatMultExpr<VT, MT>&);
  template <typename VT, typename MT> SEXP wrap(const blaze::TDVecSMatMultExpr<VT, MT>&);
  template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TSMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatDMatSchurExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatDMatSubExpr<MT1, MT2>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::TSMatDVecMultExpr<MT, VT>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatKronExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatMultExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatSchurExpr<MT1, MT2>&);
  template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatSubExpr<MT1, MT2>&);
  template <typename MT, typename VT> SEXP wrap(const blaze::TSMatSVecMultExpr<MT, VT>&);
  template <typename VT, typename MT> SEXP wrap(const blaze::TSVecDMatMultExpr<VT, MT>&);
  template <typename VT, typename MT> SEXP wrap(const blaze::TSVecSMatMultExpr<VT, MT>&);

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
