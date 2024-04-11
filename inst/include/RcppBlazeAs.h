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
    template<typename Type, bool TF>
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
    template<typename Type, bool TF>
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
    template<typename Type, bool TF>
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
    template<typename Type, bool TF>
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

    // ---------------------------- Sparse Vector Exporter ----------------------------

    /*
    // Provides only blaze::CompressedVector<Type,blaze::rowVector> export
    template<typename Type>
    class Exporter< blaze::CompressedVector<Type, blaze::rowVector> > {
      Rcpp::S4 vec;

    public:
      Exporter(SEXP x) : vec(x) {
        Rcpp::IntegerVector dims = vec.slot("Dim");
        if (dims[0] != 1) {
          throw std::invalid_argument("Not a sparse row vector.");
        }
      }

      blaze::CompressedVector<Type,blaze::rowVector> get() {
        const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;

        Rcpp::IntegerVector dims = vec.slot("Dim");
        size_t ncol = (size_t) dims[1];

        // inistialize sparse rowvector
        blaze::CompressedVector<Type, blaze::rowVector> result(ncol);

        // get the type of sparse matrix
        const char* spMatClass = vec.slot("class");

        // get data of sparse matrix
        Rcpp::Vector<RTYPE> x = vec.slot("x");

        if (strcmp(spMatClass + 2, "CMatrix") == 0) {
          Rcpp::IntegerVector p = vec.slot("p");

          size_t j = 0UL;
          for (size_t k=1UL; k < (size_t) p.size(); ++k) {
            if ((p[k] > 0UL) && (j < p[k])) {
              result(0UL, j-1UL) = x[j];
              ++j;
            }
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(vec.slot("diag"));
            if (diag == "U") {
              result(0UL, 0UL) = 1.0;
            }
          }
        } else if (strcmp(spMatClass + 2, "TMatrix") == 0) {
          Rcpp::IntegerVector j = vec.slot("j");
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            result(0UL, j[k]) = x[k];
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(vec.slot("diag"));
            if (diag == "U") {
              result(0UL, 0UL) = 1.0;
            }
          }
        } else if (strcmp(spMatClass + 2, "RMatrix") == 0) {
          Rcpp::IntegerVector j = vec.slot("j");
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            result(0UL, j[k]) = x[k];
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(vec.slot("diag"));
            if (diag == "U") {
              result(0UL, 0UL) = 1.0;
            }
          }
        } else {
          std::string spMatClassStr(spMatClass);
          Rcpp::stop(spMatClassStr + " is not supported.");
        }
        return result;
      }
    };

    // Provides only blaze::CompressedVector<Type,blaze::columnVector> export
    template<typename Type>
    class Exporter< blaze::CompressedVector<Type, blaze::columnVector> > {
      Rcpp::S4 vec;

    public:
      Exporter(SEXP x) : vec(x) {
        Rcpp::IntegerVector dims = vec.slot("Dim");
        if (dims[1] != 1) {
          throw std::invalid_argument("Not a sparse column vector.");
        }
      }

      blaze::CompressedVector<Type,blaze::columnVector> get() {
        const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;

        Rcpp::IntegerVector dims = vec.slot("Dim");
        size_t nrow = (size_t) dims[0];

        // inistialize sparse column vector
        blaze::CompressedVector<Type, blaze::columnVector> result(nrow);

        // get the type of sparse matrix
        const char* spMatClass = vec.slot("class");

        // get data of sparse matrix
        Rcpp::Vector<RTYPE> x = vec.slot("x");

        if (strcmp(spMatClass + 2, "CMatrix") == 0) {
          Rcpp::IntegerVector i = vec.slot("i");

          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            result(i[k], 0UL) = x[k];
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(vec.slot("diag"));
            if (diag == "U") {
              result(0UL, 0UL) = 1.0;
            }
          }
        } else if (strcmp(spMatClass + 2, "TMatrix") == 0) {
          Rcpp::IntegerVector i = vec.slot("i");
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            result(i[k], 0UL) = x[k];
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(vec.slot("diag"));
            if (diag == "U") {
              result(0UL, 0UL) = 1.0;
            }
          }
        } else if (strcmp(spMatClass + 2, "RMatrix") == 0) {
          Rcpp::IntegerVector j = vec.slot("j");
          Rcpp::IntegerVector p = vec.slot("p");

          size_t i = 0UL;
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            // TODO: logical problem to fix
            if (k == (size_t) p[i+1UL]) {
              ++i;
            }
            result(i, j[k]) = x[k];
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(vec.slot("diag"));
            if (diag == "U") {
              result(0UL, 0UL) = 1.0;
            }
          }
        } else {
          std::string spMatClassStr(spMatClass);
          Rcpp::stop(spMatClassStr + " is not supported.");
        }
        return result;
      }
    };
*/


  // ---------------------------- Sparse Matrix Exporter ----------------------------
  /*
    // Provides only blaze::CompressedMatrix<Type,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CompressedMatrix<Type, SO> > {
      Rcpp::S4 mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CompressedMatrix<Type, SO> get() {
        const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;

        Rcpp::IntegerVector dims = mat.slot("Dim");
        size_t nrow = (size_t) dims[0];
        size_t ncol = (size_t) dims[1];

        // initialize blaze::CompressedMatrix
        blaze::CompressedMatrix<Type,SO> result(nrow, ncol);

        // get the type of sparse matrix
        const char* spMatClass = mat.slot("class");

        // get data of sparse matrix
        Rcpp::Vector<RTYPE> x = mat.slot("x");

        if (strcmp(spMatClass + 2, "CMatrix") == 0) {
          Rcpp::IntegerVector i = mat.slot("i");
          Rcpp::IntegerVector p = mat.slot("p");

          size_t j = 0UL;
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            // TODO: logical problem to fix
            if (k == (size_t) p[j+1]) {
              ++j;
            }
            result(i[k], j) = x[k];
            if (spMatClass[1] == 's') {
              result(j, i[k]) = x[k];
            }
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(mat.slot("diag"));
            if (diag == "U") {
              auto diag(band(result, 0L));
              diag = 1.0;
            }
          }
        } else if (strcmp(spMatClass + 2, "TMatrix") == 0) {
          Rcpp::IntegerVector j = mat.slot("j");
          Rcpp::IntegerVector p = mat.slot("p");

          size_t i = 0UL;
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            // TODO: logical problem to fix
            if (k == (size_t) p[i+1]) {
              ++i;
            }
            result(i, j[k]) = x[k];
            if (spMatClass[1] == 's') {
              result(j[k], i) = x[k];
            }
          }
        } else if (strcmp(spMatClass + 2, "RMatrix") == 0) {
          Rcpp::IntegerVector i = mat.slot("i");
          Rcpp::IntegerVector j = mat.slot("j");
          for (size_t k=0UL; k < (size_t) x.size(); ++k) {
            result(i[k], j[k]) = x[k];
            if (spMatClass[1] == 's') {
              result(j[k], i[k]) = x[k];
            }
          }

          if (spMatClass[1] == 't') {
            std::string diag = Rcpp::as<std::string>(mat.slot("diag"));
            if (diag == "U") {
              auto diag(band(result, 0L));
              diag = 1.0;
            }
          }
        }  else if (strcmp(spMatClass, "ddiMatrix") == 0) {
          Rcpp::IntegerVector i = mat.slot("i");
          Rcpp::IntegerVector p = mat.slot("p");
          std::string diag = Rcpp::as<std::string>(mat.slot("diag"));

          for (size_t k=0UL; k < ncol; ++k) {
            result(k, k) = (diag == "U") ? 1.0 : x[k];
          }
        } else {
          std::string spMatClassStr(spMatClass);
          Rcpp::stop(spMatClassStr + " is not supported.");
        }
        return result;
      }
    };
    */
  } // end traits
}

#endif
