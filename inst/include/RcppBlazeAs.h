// Copyright (C) 2010 - 2024 Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C) 2017 - 2024 Ching-Chuan Chen
//
// This file is based on Rcpp.
// This file is part of RcppBlaze.
//
// RcppBlaze is free software: you can redistribute it and/or modify it
// under the terms of the 3-Clause BSD License. You should have received
// a copy of 3-Clause BSD License along with RcppBlaze.
// If not, see https://opensource.org/license/BSD-3-Clause.
#ifndef RcppBlaze__RcppBlazeAs__h
#define RcppBlaze__RcppBlazeAs__h

#include "RcppBlazeForward.h"
#include <Rcpp.h>

// most code are modified from the Exporter class definition in Rcpp / RcppArmadillo / RcppEigen
namespace Rcpp {

  namespace traits {

    // ---------------------------- Dense Vector Exporter ----------------------------
    // exporter for DynamicVector/HybridVector
    template <typename T, typename value_type>
    class BlazeVectorExporter {
    public:
      typedef value_type r_export_type;

      BlazeVectorExporter(SEXP r_obj): object(r_obj){}
      ~BlazeVectorExporter(){}

      T get() {
        R_xlen_t n = ::Rf_xlength(object);
        T result(n);
        value_type* data = blaze::data(result);
        Rcpp::internal::export_indexing<value_type*, value_type>(object, data);
        return result;
      }

    private:
      SEXP object;
    };

    template <typename Type, bool TF>
    class Exporter< blaze::DynamicVector<Type, TF> > : public BlazeVectorExporter<blaze::DynamicVector<Type, TF>, Type> {
    public:
      Exporter(SEXP x) : BlazeVectorExporter<blaze::DynamicVector<Type, TF>, Type>(x){}
    };

    template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::HybridVector<Type, N, TF, AF, PF> > : public BlazeVectorExporter<blaze::HybridVector<Type, N, TF, AF, PF>, Type> {
    public:
      Exporter(SEXP x) : BlazeVectorExporter<blaze::HybridVector<Type, N, TF, AF, PF>, Type>(x){}
    };

#define RCPPBLAZE_GET_TYPEMAP                                             \
  const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;         \
  typedef typename Rcpp::traits::storage_type<RTYPE>::type r_object_type; \
  Shield<SEXP> x(Rcpp::r_cast<RTYPE>(object));                            \
  r_object_type* y = Rcpp::internal::r_vector_start<RTYPE>(x);

#define RCPPBLAZE_VEC_COPY(__SIZE__)                                      \
  for (size_t i=0UL; i<__SIZE__; ++i) {                                   \
    result[i] = Rcpp::internal::caster<r_object_type, Type>(y[i]);        \
  }


    // exporter for StaticVector since its constructor does not need size
    template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticVector<Type, N, TF, AF, PF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::StaticVector<Type, N, TF, AF, PF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        size_t n = (size_t) Rf_xlength(x);
        if (n != N) {
          throw std::invalid_argument("Dimension of vector is not match.");
        }

        blaze::StaticVector<Type, N, TF, AF, PF> result;
        RCPPBLAZE_VEC_COPY(n);
        return result;
      }

    private:
      SEXP object;
    };

    // Provides only blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> export
    template <typename Type, bool TF>
    class Exporter< blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        size_t n = (size_t) Rf_xlength(x);
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[n]);
        blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> result(data.get(), n);
        RCPPBLAZE_VEC_COPY(n);
        return result;
      }

    private:
      SEXP object;
    };

    // Provides only blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> export
    template <typename Type, bool TF>
    class Exporter< blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        size_t n = (size_t) Rf_xlength(x);
        size_t paddedSize = blaze::nextMultiple<size_t>(n, blaze::SIMDTrait<Type>::size);
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[paddedSize]);
        blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> result(data.get(), n, paddedSize);
        RCPPBLAZE_VEC_COPY(n);
        return result;
      }

    private:
      SEXP object;
    };

    // Provides only blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> export (not supported)
    template <typename Type, bool TF>
    class Exporter< blaze::CustomVector<Type, blaze::aligned, blaze::unpadded,TF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        size_t n = (size_t) Rf_xlength(x);
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(n));
        blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> result(data.get(), n);
        RCPPBLAZE_VEC_COPY(n);
        return result;
      }

    private:
      SEXP object;
    };

    // Provides only blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> export (not supported)
    template <typename Type, bool TF>
    class Exporter< blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        size_t n = (size_t) Rf_xlength(x);
        size_t paddedSize = blaze::nextMultiple<size_t>(n, blaze::SIMDTrait<Type>::size);
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(paddedSize));
        blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> result(data.get(), n, paddedSize);
        RCPPBLAZE_VEC_COPY(n);
        return result;
      }

    private:
      SEXP object;
    };

#undef RCPPBLAZE_VEC_COPY

    // ---------------------------- Dense Matrix Exporter ----------------------------
#define RCPPBLAZE_MATRIX_COPY(__ROWS__, __COLS__)                                       \
  for (size_t j=0UL; j<__COLS__; ++j) {                                                 \
    for (size_t i=0UL; i<__ROWS__; ++i) {                                               \
      result(i, j) = Rcpp::internal::caster<r_object_type, Type>(y[j*__ROWS__+i]); \
    }}

#define RCPPBLAZE_GET_MATRIX_DIMS(__ROWS__, __COLS__)                    \
  Shield<SEXP> obj_dims(::Rf_getAttrib(object, R_DimSymbol));            \
  if (Rf_isNull(obj_dims) || ::Rf_length(obj_dims) != 2) {               \
    throw ::Rcpp::not_a_matrix();                                        \
  }                                                                      \
  int* dims = INTEGER(obj_dims);                                         \
  size_t __ROWS__ = (size_t) dims[0], __COLS__ = (size_t) dims[1];

#define RCPPBLAZE_DIM_MISMATCH_EXCEPTION(__ROWS__, __CLS_ROWS__, __COLS__, __CLS_COLS__)  \
  if ((__ROWS__ != __CLS_ROWS__) || (__COLS__ != __CLS_COLS__)) {                         \
    throw std::invalid_argument("Dimension of matrix is not match.");                     \
  }

    // exporter for DynamicMatrix
    template <typename Type, bool SO>
    class Exporter< blaze::DynamicMatrix<Type, SO> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::DynamicMatrix<Type, SO> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        blaze::DynamicMatrix<Type, SO> result(m, n);
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

    // exporter for HybridMatrix
    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::HybridMatrix<Type, M, N, SO, AF, PF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::HybridMatrix<Type, M, N, SO, AF, PF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        RCPPBLAZE_DIM_MISMATCH_EXCEPTION(m, M, n, N);
        blaze::HybridMatrix<Type, M, N, SO, AF, PF> result(m, n);
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

    // exporter for StaticMatrix since its constructor does not need size
    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticMatrix<Type, M, N, SO, AF, PF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::StaticMatrix<Type, M, N, SO, AF, PF> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        RCPPBLAZE_DIM_MISMATCH_EXCEPTION(m, M, n, N);
        blaze::StaticMatrix<Type, M, N, SO, AF, PF> result;
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

    // exporter for blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO>
    template <typename Type, bool SO>
    class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[m*n]);
        blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> result(data.get(), m, n);
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

    // exporter for blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO>
    template <typename Type, bool SO>
    class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        size_t simdTypeSize = blaze::SIMDTrait<Type>::size, paddedSize, matSize;
        if (SO == blaze::rowMajor) {
          paddedSize = blaze::nextMultiple<size_t>(n, simdTypeSize);
          matSize = m * paddedSize;
        } else {
          paddedSize = blaze::nextMultiple<size_t>(m, simdTypeSize);
          matSize = paddedSize * n;
        }
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[matSize]);
        blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO> result(data.get(), m, n, paddedSize);
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

    // exporter for blaze::CustomMatrix<Type, blaze::aligned, blaze::unpadded, SO>
    template <typename Type, bool SO>
    class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::unpadded, SO> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj): object(r_obj){}
      ~Exporter(){}

      blaze::CustomMatrix<Type, blaze::aligned, blaze::unpadded, SO> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        size_t simdTypeSize = blaze::SIMDTrait<Type>::size, paddedSize, matSize;
        if (SO == blaze::rowMajor) {
          paddedSize = blaze::nextMultiple<size_t>(n, simdTypeSize);
          matSize = m * paddedSize;
        } else {
          paddedSize = blaze::nextMultiple<size_t>(m, simdTypeSize);
          matSize = paddedSize * n;
        }
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(matSize));
        blaze::CustomMatrix<Type, blaze::aligned, blaze::unpadded, SO> result(data.get(), m, n, paddedSize);
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

    // exporter for blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO>
    template <typename Type, bool SO>
    class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO> get() {
        RCPPBLAZE_GET_TYPEMAP;
        RCPPBLAZE_GET_MATRIX_DIMS(m, n);
        size_t simdTypeSize = blaze::SIMDTrait<Type>::size, paddedSize, matSize;
        if (SO == blaze::rowMajor) {
          paddedSize = blaze::nextMultiple<size_t>(n, simdTypeSize);
          matSize = m * paddedSize;
        } else {
          paddedSize = blaze::nextMultiple<size_t>(m, simdTypeSize);
          matSize = paddedSize * n;
        }
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(matSize));
        blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO> result(data.get(), m, n, paddedSize);
        RCPPBLAZE_MATRIX_COPY(m, n);
        return result;
      }

    private:
      SEXP object;
    };

#undef RCPPBLAZE_GET_TYPEMAP
#undef RCPPBLAZE_GET_MATRIX_DIMS
#undef RCPPBLAZE_MATRIX_COPY
#undef RCPPBLAZE_DIM_MISMATCH_EXCEPTION

    // ---------------------------- Sparse Vector Exporter ----------------------------

#define RCPPBLAZE_GET_S4_OBJ_MAP                                                     \
    const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;                  \
    if ((RTYPE != INTSXP) && (RTYPE != REALSXP)) {                                   \
      Rcpp::stop("SparseVector/SparseMatrix only supports int, float and double!");  \
    }                                                                                \
    Rcpp::IntegerVector dims = object.slot("Dim");                                   \
    int m = dims[0], n = dims[1];

#define RCPPBLAZE_S4_MAT_CSC                                     \
  Rcpp::IntegerVector i = object.slot("i"), p = object.slot("p");

#define RCPPBLAZE_S4_MAT_CSR                                     \
  Rcpp::IntegerVector j = object.slot("j"), p = object.slot("p");

#define RCPPBLAZE_S4_MAT_COO                                     \
  Rcpp::IntegerVector i = object.slot("i"), j = object.slot("j");

    // Provides only blaze::CompressedVector<Type, TF> export
    template <typename Type, bool TF>
    class Exporter< blaze::CompressedVector<Type, TF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::CompressedVector<Type, TF> get() {
        RCPPBLAZE_GET_S4_OBJ_MAP;
        size_t size;
        if (TF == blaze::rowVector) {
          if (m == 1) {
            size = (size_t) n;
          } else {
            Rcpp::stop("The number of rows should be 1 for CompressedVector<Type, rowVector>!");
          }
        } else {
          if (n == 1) {
            size = (size_t) m;
          } else {
            Rcpp::stop("The number of columns should be 1 for CompressedVector<Type, rowVector>!");
          }
        }

        std::string matrixClass = Rcpp::as<std::string>(object.slot("class"));
        Rcpp::Vector<RTYPE> x = object.slot("x");
        size_t nonZeroSize = (size_t) x.size();
        blaze::CompressedVector<Type, TF> result(size, nonZeroSize);
        size_t u, v;
        int cnt;
        if (matrixClass == "dgCMatrix" || object.is("dgCMatrix")) {
          RCPPBLAZE_S4_MAT_CSC;
          size_t pSize = (size_t) p.size();
          v = 0UL;
          for (u=1UL; u < pSize; ++u) {
            if (p[u] > p[u-1]) {
              cnt = 0;
              while (p[u] > p[u-1] + cnt) {
                if (TF == blaze::rowVector) {
                  result.append(u-1, x[v]);
                } else {
                  result.append((size_t)i[v], x[v]);
                }
                ++v;
                ++cnt;
              }
            }
          }
        } else if (matrixClass == "dgTMatrix" || object.is("dgTMatrix")) {
          RCPPBLAZE_S4_MAT_COO;
          for (u=0UL; u < nonZeroSize; ++u) {
            if (TF == blaze::rowVector) {
              result.append((size_t)j[v], x[v]);
            } else {
              result.append((size_t)i[v], x[v]);
            }
          }
        } else if (matrixClass == "dgRMatrix" || object.is("dgRMatrix")) {
          RCPPBLAZE_S4_MAT_CSR;
          size_t pSize = (size_t) p.size();
          v = 0UL;
          for (u=1UL; u < pSize; ++u) {
            if (p[u] > p[u-1]) {
              cnt = 0;
              while (p[u] > p[u-1] + cnt) {
                if (TF == blaze::rowVector) {
                  result.append((size_t)j[v], x[v]);
                } else {
                  result.append(u-1, x[v]);
                }
                ++v;
                ++cnt;
              }
            }
          }
        } else {
          Rcpp::stop(matrixClass + " is not supported.");
        }
        return result;
      }

    private:
      Rcpp::S4 object;
    };

    template <typename Type, bool TF>
    class Exporter< blaze::ZeroVector<Type, TF> > : public Exporter< blaze::CompressedVector<Type, TF> > {
    public:
      Exporter(SEXP x) : Exporter< blaze::CompressedVector<Type, TF> >(x){}
    };

  // ---------------------------- Sparse Matrix Exporter ----------------------------
#define RCPPBLAZE_COLUMN_MAJOR_SPARSE_MATRIX_COPY                \
    for (size_t u=0UL; u<vi_pair_vec.size(); ++u) {              \
      result.reserve(u, vi_pair_vec[u].size());                  \
      for (auto &el: vi_pair_vec[u]) {                           \
        result.append(std::get<0>(el), u, std::get<1>(el));      \
      }                                                          \
      result.finalize(u);                                        \
    }

#define RCPPBLAZE_ROW_MAJOR_SPARSE_MATRIX_COPY                 \
  for (size_t u=0UL; u<vi_pair_vec.size(); ++u) {              \
    result.reserve(u, vi_pair_vec[u].size());                  \
    for (auto &el: vi_pair_vec[u]) {                           \
      result.append(u, std::get<0>(el), std::get<1>(el));      \
    }                                                          \
    result.finalize(u);                                        \
  }

#define RCPPBLAZE_MATRIX_COPY                                  \
  if (SO == blaze::rowMajor) {                                 \
    RCPPBLAZE_ROW_MAJOR_SPARSE_MATRIX_COPY;                    \
  } else {                                                     \
    RCPPBLAZE_COLUMN_MAJOR_SPARSE_MATRIX_COPY;                 \
  }

#define RCPPBLAZE_S4_MAT_DIAG                                    \
  std::string diag = Rcpp::as<std::string>(object.slot("diag")); \
  if (diag == "U") {                                             \
    for (size_t i=0UL; i<(size_t)std::min(m, n); ++i) {          \
      result(i, i) = 1.0;                                        \
    }}

#define RCPPBLAZE_SYM_COPY(__START__, __END__)                   \
  for (size_t i=0UL; i<(size_t)m; ++i) {                         \
    for (size_t j=__START__; j < (size_t)__END__; ++j) {         \
      result(j, i) = result(i, j);                               \
    }}

#define RCPPBLAZE_S4_MAT_UPLO                                    \
  std::string uplo = Rcpp::as<std::string>(object.slot("uplo")); \
  if (uplo == "U") {                                             \
    RCPPBLAZE_SYM_COPY(i, n);                                    \
  } else {                                                       \
    RCPPBLAZE_SYM_COPY(0UL, i+1);                                \
  }


    // Provides only blaze::CompressedMatrix<Type,SO> export
    template <typename Type, bool SO >
    class Exporter< blaze::CompressedMatrix<Type, SO> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      using vi_pair = std::tuple<size_t, double>;
      std::vector<std::vector<vi_pair>> getViPair(
          bool targetSO,
          bool sourceSO,
          Rcpp::S4 object
      ) {
        std::string idxSym = (sourceSO == blaze::rowMajor)?"j":"i";
        Rcpp::IntegerVector idx = object.slot(idxSym), p = object.slot("p");
        Rcpp::NumericVector x = object.slot("x");
        Rcpp::IntegerVector dims = object.slot("Dim");

        int v = 0, cnt;
        size_t pair_size = (targetSO == blaze::rowMajor)?dims[0]:dims[1];
        std::vector<std::vector<vi_pair>> vi_pair_vec(pair_size);
        for (int u=1; u < p.size(); ++u) {
          if (p[u] > p[u-1]) {
            cnt = p[u] - p[u-1];
            vi_pair_vec[u-1].reserve(cnt);
            for (int t=0; t<cnt; ++t) {
              if (sourceSO != targetSO) {
                vi_pair_vec[idx[v]].push_back(vi_pair(u-1, x[v]));
              } else {
                vi_pair_vec[u-1].push_back(vi_pair(idx[v], x[v]));
              }
              ++v;
            }
          }
        }
        return vi_pair_vec;
      }

      std::vector<std::vector<vi_pair>> getViPairForTSparse(
          bool targetSO,
          Rcpp::S4 object
      ) {
        Rcpp::IntegerVector i = object.slot("i"), j = object.slot("j");
        Rcpp::NumericVector x = object.slot("x");
        Rcpp::IntegerVector dims = object.slot("Dim");

        size_t pair_size = (targetSO == blaze::rowMajor)?dims[0]:dims[1];
        std::vector<std::vector<vi_pair>> vi_pair_vec(pair_size);
        for (size_t u=0UL; u<(size_t)i.size(); ++u) {
          if (targetSO == blaze::rowMajor) {
            vi_pair_vec[i[u]].push_back(vi_pair(j[u], x[u]));
          } else {
            vi_pair_vec[j[u]].push_back(vi_pair(i[u], x[u]));
          }
        }
        return vi_pair_vec;
      }

      blaze::CompressedMatrix<Type, SO> get() {
        RCPPBLAZE_GET_S4_OBJ_MAP;

        Rcpp::Vector<RTYPE> x = object.slot("x");
        blaze::CompressedMatrix<Type, SO> result(m, n);
        std::string matrixClass = Rcpp::as<std::string>(object.slot("class"));

        if (matrixClass == "dgCMatrix" || object.is("dgCMatrix")) {
          RCPPBLAZE_S4_MAT_CSC;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPair(SO, blaze::columnMajor, object);
          RCPPBLAZE_MATRIX_COPY;
        } else if (matrixClass == "dsCMatrix" || object.is("dsCMatrix")) {
          RCPPBLAZE_S4_MAT_CSC;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPair(SO, blaze::columnMajor, object);
          RCPPBLAZE_MATRIX_COPY;
          RCPPBLAZE_S4_MAT_UPLO;
        } else if (matrixClass == "dtCMatrix" || object.is("dtCMatrix")) {
          RCPPBLAZE_S4_MAT_CSC;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPair(SO, blaze::columnMajor, object);
          RCPPBLAZE_MATRIX_COPY;
          RCPPBLAZE_S4_MAT_DIAG;
        } else if (matrixClass == "dgTMatrix" || object.is("dgTMatrix")) {
          RCPPBLAZE_S4_MAT_COO;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPairForTSparse(SO, object);
          RCPPBLAZE_MATRIX_COPY;
        } else if (matrixClass == "dsTMatrix" || object.is("dsTMatrix")) {
          RCPPBLAZE_S4_MAT_COO;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPairForTSparse(SO, object);
          RCPPBLAZE_MATRIX_COPY;
          RCPPBLAZE_S4_MAT_UPLO;
        } else if (matrixClass == "dtTMatrix" || object.is("dtTMatrix")) {
          RCPPBLAZE_S4_MAT_COO;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPairForTSparse(SO, object);
          RCPPBLAZE_MATRIX_COPY;
          RCPPBLAZE_S4_MAT_DIAG;
        } else if (matrixClass == "dgRMatrix" || object.is("dgRMatrix")) {
          RCPPBLAZE_S4_MAT_CSR;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPair(SO, blaze::rowMajor, object);
          RCPPBLAZE_MATRIX_COPY;
        } else if (matrixClass == "dsRMatrix" || object.is("dsRMatrix")) {
          RCPPBLAZE_S4_MAT_CSR;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPair(SO, blaze::rowMajor, object);
          RCPPBLAZE_MATRIX_COPY;
          RCPPBLAZE_S4_MAT_UPLO;
        } else if (matrixClass == "dtRMatrix" || object.is("dtRMatrix")) {
          RCPPBLAZE_S4_MAT_CSR;
          std::vector<std::vector<vi_pair>> vi_pair_vec = getViPair(SO, blaze::rowMajor, object);
          RCPPBLAZE_MATRIX_COPY;
          RCPPBLAZE_S4_MAT_DIAG;
        } else if (matrixClass == "ddiMatrix" || object.is("ddiMatrix")) {
          Rcpp::IntegerVector dims = object.slot("Dim");
          std::string diag = Rcpp::as<std::string>(object.slot("diag"));
          if (diag == "U") {
            for (size_t i=0UL; i<(size_t)dims[0]; ++i) {
              result.reserve(i, 1);
              result.append(i, i, 1.0);
              result.finalize(i);
            }
          } else {
            Rcpp::NumericVector x = object.slot("x");
            for (size_t i=0UL; i<(size_t)dims[0]; ++i) {
              result.reserve(i, 1);
              result.append(i, i, x[i]);
              result.finalize(i);
            }
          }
        } else {
          Rcpp::stop(matrixClass + " is not supported.");
        }
        return result;
      }

    private:
      Rcpp::S4 object;
    };

    template <typename Type, bool SO>
    class Exporter< blaze::ZeroMatrix<Type, SO> > : public Exporter< blaze::CompressedMatrix<Type, SO> > {
    public:
      Exporter(SEXP x) : Exporter< blaze::CompressedMatrix<Type, SO> >(x){}
    };

#undef RCPPBLAZE_GET_S4_OBJ_MAP
#undef RCPPBLAZE_S4_MAT_CSC
#undef RCPPBLAZE_S4_MAT_CSR
#undef RCPPBLAZE_S4_MAT_COO
#undef RCPPBLAZE_COLUMN_MAJOR_SPARSE_MATRIX_COPY
#undef RCPPBLAZE_ROW_MAJOR_SPARSE_MATRIX_COPY
#undef RCPPBLAZE_MATRIX_COPY
#undef RCPPBLAZE_S4_MAT_DIAG
#undef RCPPBLAZE_S4_MAT_UPLO
#undef RCPPBLAZE_SYM_COPY
  } // end traits
}

#endif
