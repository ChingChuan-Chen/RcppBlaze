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

      BlazeVectorExporter(SEXP r_obj) : object(r_obj){}
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

#define RCPPBLAZE_GET_TYPEMAP_SIZE                                        \
  const int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;         \
  typedef typename Rcpp::traits::storage_type<RTYPE>::type r_object_type; \
  Shield<SEXP> x(Rcpp::r_cast<RTYPE>(object));                            \
  r_object_type* y = Rcpp::internal::r_vector_start<RTYPE>(x);            \
  size_t size = (size_t) Rf_xlength(x);

#define RCPPBLAZE_VEC_COPY                                                \
  for (size_t i=0UL; i<size; ++i) {                                       \
    result[i] = Rcpp::internal::caster<r_object_type, Type>(y[i]);        \
  }


    // exporter for StaticVector since its constructor does not need size
    template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticVector<Type, N, TF, AF, PF> > {
    public:
      typedef Type r_export_type;

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::StaticVector<Type, N, TF, AF, PF> get() {
        RCPPBLAZE_GET_TYPEMAP_SIZE;
        if (size != N) {
          throw std::invalid_argument("Dimension of vector is not match.");
        }

        blaze::StaticVector<Type, N, TF, AF, PF> result;
        RCPPBLAZE_VEC_COPY;
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

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP_SIZE;
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[size]);
        blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> result(data.get(), size);
        RCPPBLAZE_VEC_COPY;
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

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP_SIZE;
        size_t paddedSize = blaze::nextMultiple<size_t>(size, blaze::SIMDTrait<Type>::size);
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[paddedSize]);
        blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> result(data.get(), size, paddedSize);
        RCPPBLAZE_VEC_COPY;
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

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP_SIZE;
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(size));
        blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> result(data.get(), size);
        RCPPBLAZE_VEC_COPY;
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

      Exporter(SEXP r_obj) : object(r_obj){}
      ~Exporter(){}

      blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> get() {
        RCPPBLAZE_GET_TYPEMAP_SIZE;
        size_t paddedSize = blaze::nextMultiple<size_t>(size, blaze::SIMDTrait<Type>::size);
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(paddedSize));
        blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> result(data.get(), size, paddedSize);
        RCPPBLAZE_VEC_COPY;
        return result;
      }

    private:
      SEXP object;
    };

#undef RCPPBLAZE_VEC_COPY

    // ---------------------------- Dense Matrix Exporter ----------------------------
#define RCPPBLAZE_MATRIX_COPY(__IDX1__, __IDX2__, __END1__, __END2__)    \
    for (size_t __IDX1__=0UL; __IDX1__ < __END1__; ++__IDX1__) {           \
      for (size_t __IDX2__=0UL; __IDX2__ < __END2__; ++__IDX2__) {         \
        result(i, j) = Rcpp::internal::caster<value_t, Type>(mat(i, j)); \
      }}

    // exporter for DynamicMatrix
    template <typename Type, bool SO>
    class Exporter< blaze::DynamicMatrix<Type, SO> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : mat(x) {}
      blaze::DynamicMatrix<Type, SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;

        size_t m = mat.nrow(), n = mat.ncol();
        blaze::DynamicMatrix<Type, SO> result(m, n);
        if (SO == blaze::rowMajor) {
          RCPPBLAZE_MATRIX_COPY(i, j, m, n);
        } else {
          RCPPBLAZE_MATRIX_COPY(j, i, n, m);
        }
        return result;
      }
    };

    // exporter for HybridMatrix
    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::HybridMatrix<Type, M, N, SO, AF, PF> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : mat(x) {}
      blaze::HybridMatrix<Type, M, N, SO, AF, PF> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;

        size_t m = mat.nrow(), n = mat.ncol();
        if ((m != M) || (n != N)) {
          throw std::invalid_argument("Dimension of matrix is not match.");
        }

        blaze::HybridMatrix<Type, M, N, SO, AF, PF> result(m, n);
        if (SO == blaze::rowMajor) {
          RCPPBLAZE_MATRIX_COPY(i, j, m, n);
        } else {
          RCPPBLAZE_MATRIX_COPY(j, i, n, m);
        }
        return result;
      }
    };

    // exporter for StaticMatrix since its constructor does not need size
    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticMatrix<Type, M, N, SO, AF, PF> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : mat(x) {}
      blaze::StaticMatrix<Type, M, N, SO, AF, PF> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;

        size_t m = mat.nrow(), n = mat.ncol();
        if ((m != M) || (n != N)) {
          throw std::invalid_argument("Dimension of matrix is not match.");
        }

        blaze::StaticMatrix<Type, M, N, SO, AF, PF> result;
        if (SO == blaze::rowMajor) {
          RCPPBLAZE_MATRIX_COPY(i, j, m, n);
        } else {
          RCPPBLAZE_MATRIX_COPY(j, i, n, m);
        }
        return result;
      }
    };

    // exporter for blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO>
    template <typename Type, bool SO>
    class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        size_t m = mat.nrow(), n = mat.ncol();
        size_t size = m*n;
        std::unique_ptr<Type[], blaze::ArrayDelete> data(new Type[size]);

        blaze::CustomMatrix<Type, blaze::unaligned, blaze::unpadded, SO> result(data.get(), m, n);
        if (SO == blaze::rowMajor) {
          RCPPBLAZE_MATRIX_COPY(i, j, m, n);
        } else {
          RCPPBLAZE_MATRIX_COPY(j, i, n, m);
        }
        return result;
      }
    };

#undef RCPPBLAZE_MATRIX_COPY
    /*
    // Provides only blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded, SO> get() {
        blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded, SO> result(
            reinterpret_cast<Type*>(mat.begin()),
            (size_t) mat.nrow(),
            (size_t) mat.ncol()
        );
        return result ;
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::unpadded, SO> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        if (SO == blaze::rowMajor) {
          size_t actualCols = (size_t) mat.ncol() + Type::size - (size_t) mat.ncol() % Type::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> result(
              blaze::allocate<Type>( actualCols * (size_t)mat.nrow() ), (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualCols, blaze::Deallocate() );
          for (size_t i=0UL; i<result.rows(); ++i) {
            for (size_t j=0UL; j<result.columns(); ++j) {
              result(i ,j) = Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
          return result;
        } else {
          size_t actualRows = (size_t) mat.nrow() + Type::size - (size_t) mat.nrow() % Type::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::unpadded,SO> result(
              blaze::allocate<Type>( actualRows * (size_t)mat.ncol() ), (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualRows, blaze::Deallocate() );
          for (size_t j=0UL; j<result.columns(); ++j) {
            for (size_t i=0UL; i<result.rows(); ++i) {
              result(i ,j) = Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
          return result;
        }
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        if (SO == blaze::rowMajor) {
          size_t actualCols = (size_t) mat.ncol() + Type::size - (size_t) mat.ncol() % Type::size;
          blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> result(
              new Type[actualCols * (size_t)mat.nrow()],
                      (size_t) mat.nrow(),
                      (size_t) mat.ncol(),
                      actualCols,
                      blaze::ArrayDelete()
          );
          for (size_t i=0UL; i<result.rows(); ++i) {
            for (size_t j=0UL; j<result.columns(); ++j) {
              result(i ,j) = Rcpp::internal::caster<value_t, Type>(mat(i ,j));
            }
          }
          return result;
        } else {
          size_t actualRows = (size_t) mat.nrow() + Type::size - (size_t) mat.nrow() % Type::size;
          blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> result(
              new Type[actualRows * (size_t)mat.ncol()],
                      (size_t) mat.nrow(),
                      (size_t) mat.ncol(),
                      actualRows,
                      blaze::ArrayDelete()
          );
          for (size_t j=0UL; j<result.columns(); ++j) {
            for (size_t i=0UL; i<result.rows(); ++i) {
              result(i ,j) = Rcpp::internal::caster<value_t, Type>(mat(i ,j));
            }
          }
          return result;
        }
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO> > {
      const static int RTYPE = Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        if (SO == blaze::rowMajor) {
          size_t actualCols = (size_t) mat.ncol() + Type::size - (size_t) mat.ncol() % Type::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> result(
              blaze::allocate<Type>( actualCols * (size_t)mat.nrow() ),
              (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualCols, blaze::Deallocate()
          );
          for (size_t i=0UL; i<result.rows(); ++i) {
            for (size_t j=0UL; j<result.columns(); ++j) {
              result(i ,j) = Rcpp::internal::caster<value_t, Type>(mat(i ,j));
            }
          }
          return result;
        } else {
          size_t actualRows = (size_t) mat.nrow() + Type::size - (size_t) mat.nrow() % Type::size;
          blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> result(
              blaze::allocate<Type>(actualRows * (size_t)mat.ncol()),
              (size_t) mat.nrow(),
              (size_t) mat.ncol(), actualRows, blaze::Deallocate()
          );
          for (size_t j=0UL; j<result.columns(); ++j) {
            for (size_t i=0UL; i<result.rows(); ++i) {
              result(i ,j) = Rcpp::internal::caster<value_t, Type>(mat(i ,j));
            }
          }
          return result;
        }
      }
    };
     */


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
