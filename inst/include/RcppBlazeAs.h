// Copyright (C)  2013 - 2021  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2017 - 2021  Serguei Sokol
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
// Copyright (C)  2017 - 2024  Ching-Chuan Chen
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

#include "RcppBlazeForward.h"
#include <Rcpp.h>

// most code are modified from the Exporter class definition in Rcpp / RcppArmadillo / RcppEigen
namespace Rcpp {

  namespace traits {

    // exporter for DynamicVector/HybridVector
    template <typename T, typename value_type>
    class BlazeVectorExporter {
    public:
      typedef value_type r_export_type;

      BlazeVectorExporter(SEXP x) : object(x){}
      ~BlazeVectorExporter(){}

      T get() {
        R_xlen_t n = ::Rf_xlength(object);
        T result(n);
        value_type* data = blaze::data(result);
        ::Rcpp::internal::export_indexing<value_type*, value_type>(object, data);
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

    // exporter for StaticVector since its constructor does not need size
    template <typename T, typename value_type>
    class BlazeStaticVectorExporter {
    public:
      typedef value_type r_export_type;

      BlazeStaticVectorExporter(SEXP x) : object(x){}
      ~BlazeStaticVectorExporter(){}

      T get() {
        T result;
        value_type* data = blaze::data(result);
        ::Rcpp::internal::export_indexing<value_type*, value_type>(object, data);
        return result;
      }

    private:
      SEXP object;
    };

    template <typename Type, size_t N, bool TF, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::StaticVector<Type, N, TF, AF, PF> > : public BlazeStaticVectorExporter<blaze::StaticVector<Type, N, TF, AF, PF>, Type> {
    public:
      Exporter(SEXP x) : BlazeStaticVectorExporter<blaze::StaticVector<Type, N, TF, AF, PF>, Type>(x){}
    };

    // Provides only blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> export
    template<typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Vector<RTYPE> vec;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        size_t n = (size_t) vec.size();
        std::vector<Type> data(n);

        blaze::CustomVector<Type, blaze::unaligned, blaze::unpadded, TF> result(data.data(), n);
        for (size_t i=0UL; i < (size_t) vec.size(); ++i) {
          result[i] = ::Rcpp::internal::caster<value_t, Type>(vec[i]);
        }
        return result;
      }
    };

    // Provides only blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> export
    template<typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Vector<RTYPE> vec;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        size_t n = (size_t) vec.size();
        size_t simdTypeSize = blaze::SIMDTrait<Type>::size;
        size_t paddedSize = blaze::nextMultiple<size_t>(n, simdTypeSize);
        std::unique_ptr<Type[], blaze::Deallocate> data(new Type[paddedSize]);

        blaze::CustomVector<Type, blaze::unaligned, blaze::padded, TF> result(data.get(), n, paddedSize);
        for (size_t i=0UL; i < (size_t) vec.size(); ++i) {
          result[i] = ::Rcpp::internal::caster<value_t, Type>(vec[i]);
        }
        return result;
      }
    };

    // Provides only blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> export (not supported)
    template<typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type, blaze::aligned, blaze::unpadded,TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Vector<RTYPE> vec;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> get() {
        Rcpp::stop("Aligned CustomVector is not supported");

        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        size_t n = (size_t) vec.size();
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(n));

        blaze::CustomVector<Type, blaze::aligned, blaze::unpadded, TF> result(data.get(), n);
        for (size_t i=0UL; i < (size_t) vec.size(); ++i) {
          result[i] = ::Rcpp::internal::caster<value_t, Type>(vec[i]);
        }
        return result;
      }
    };


    // Provides only blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> export (not supported)
    template<typename Type, bool TF >
    class Exporter< blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      Rcpp::Vector<RTYPE> vec;

    public:
      typedef Type r_export_type;
      Exporter(SEXP x) : vec(x) {}
      blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> get() {
        Rcpp::stop("Aligned CustomVector is not supported");

        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        size_t n = (size_t) vec.size();
        size_t simdTypeSize = blaze::SIMDTrait<Type>::size;
        size_t paddedSize = 16; // blaze::nextMultiple<size_t>(n, simdTypeSize);
        std::unique_ptr<Type[], blaze::Deallocate> data(blaze::allocate<Type>(paddedSize));

        blaze::CustomVector<Type, blaze::aligned, blaze::padded, TF> result(data.get(), n, paddedSize);
        for (size_t i=0UL; i < (size_t) vec.size(); ++i) {
          result[i] = ::Rcpp::internal::caster<value_t, Type>(vec[i]);
        }
        return result;
      }
    };


    /*
    template <typename T, typename value_type>
    class BlazeCustomVectorExporter {
    public:
      typedef value_type r_export_type;

      BlazeCustomVectorExporter(SEXP x) : object(x){}
      ~BlazeCustomVectorExporter(){}

      T get() {
        size_t n = (size_t) ::Rf_xlength(object);
        std::vector<value_type> vec = Rcpp::as< std::vector<value_type> >(object);
        T result(&vec[0], n);
        return result;
      }

    private:
      SEXP object;
    };

    template <typename Type, blaze::AlignmentFlag AF, blaze::PaddingFlag PF, bool TF>
    class Exporter< blaze::CustomVector<Type, AF, PF, TF> > : public BlazeCustomVectorExporter<blaze::CustomVector<Type, AF, PF, TF>, Type> {
    public:
      Exporter(SEXP x) : BlazeCustomVectorExporter<blaze::CustomVector<Type, AF, PF, TF>, Type>(x){}
    };
     */



    // exporter for DynamicMatrix/HybridMatrix
    template <typename T, typename value_type>
    class MatrixExporterForBlaze {
    public:
      typedef value_type r_export_type;

      MatrixExporterForBlaze(SEXP x) : object(x){}
      ~MatrixExporterForBlaze(){}

      T get() {
        Shield<SEXP> dims(::Rf_getAttrib( object, R_DimSymbol));
        if (Rf_isNull(dims) || ::Rf_length(dims) != 2) {
          throw ::Rcpp::not_a_matrix();
        }
        int* dims_ = INTEGER(dims);
        T result(dims_[0], dims_[1]);
        value_type* data = blaze::data(result);
        ::Rcpp::internal::export_indexing<value_type*, value_type>(object, data);
        return result;
      }

    private:
      SEXP object;
    };

    template<typename Type, bool SO>
    class Exporter< blaze::DynamicMatrix<Type, SO> > : public MatrixExporterForBlaze<blaze::DynamicMatrix<Type, SO>, Type> {
    public:
      Exporter(SEXP x) : MatrixExporterForBlaze<blaze::DynamicMatrix<Type, SO>, Type>(x){}
    };

    template <typename Type, size_t M, size_t N, bool SO, blaze::AlignmentFlag AF, blaze::PaddingFlag PF>
    class Exporter< blaze::HybridMatrix<Type, M, N, SO, AF, PF> > : public MatrixExporterForBlaze<blaze::HybridMatrix<Type, M, N, SO, AF, PF>, Type> {
    public:
      Exporter(SEXP x) : MatrixExporterForBlaze<blaze::HybridMatrix<Type, M, N, SO, AF, PF>, Type>(x){}
    };

    /*


    // Provides only blaze::CompressedVector<Type,blaze::rowVector> export
    template<typename Type>
    class Exporter< blaze::CompressedVector<Type, blaze::rowVector> > {
      ::Rcpp::S4 vec;

    public:
      Exporter(SEXP x) : vec(x) {
        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        if (dims[0] != 1) {
          throw std::invalid_argument("Not a sparse row vector.");
        }
      }

      blaze::CompressedVector<Type,blaze::rowVector> get() {
        const int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;

        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        size_t ncol = (size_t) dims[1];

        // inistialize sparse rowvector
        blaze::CompressedVector<Type, blaze::rowVector> result(ncol);

        // get the type of sparse matrix
        const char* spMatClass = vec.slot("class");

        // get data of sparse matrix
        ::Rcpp::Vector<RTYPE> x = vec.slot("x");

        if (strcmp(spMatClass + 2, "CMatrix") == 0) {
          ::Rcpp::IntegerVector p = vec.slot("p");

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
          ::Rcpp::IntegerVector j = vec.slot("j");
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
          ::Rcpp::IntegerVector j = vec.slot("j");
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
      ::Rcpp::S4 vec;

    public:
      Exporter(SEXP x) : vec(x) {
        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        if (dims[1] != 1) {
          throw std::invalid_argument("Not a sparse column vector.");
        }
      }

      blaze::CompressedVector<Type,blaze::columnVector> get() {
        const int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;

        ::Rcpp::IntegerVector dims = vec.slot("Dim");
        size_t nrow = (size_t) dims[0];

        // inistialize sparse column vector
        blaze::CompressedVector<Type, blaze::columnVector> result(nrow);

        // get the type of sparse matrix
        const char* spMatClass = vec.slot("class");

        // get data of sparse matrix
        ::Rcpp::Vector<RTYPE> x = vec.slot("x");

        if (strcmp(spMatClass + 2, "CMatrix") == 0) {
          ::Rcpp::IntegerVector i = vec.slot("i");

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
          ::Rcpp::IntegerVector i = vec.slot("i");
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
          ::Rcpp::IntegerVector j = vec.slot("j");
          ::Rcpp::IntegerVector p = vec.slot("p");

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

    // Provides only blaze::DynamicMatrix<Type,SO> export
    /*
    template<typename Type, bool SO >
    class Exporter< blaze::DynamicMatrix<Type, SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::DynamicMatrix<Type,SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        blaze::DynamicMatrix<Type,SO> result( (size_t) mat.nrow(), (size_t) mat.ncol() );
        if (SO == blaze::rowMajor) {
          for (size_t i=0UL; i<result.rows(); ++i) {
            for (size_t j=0UL; j<result.columns(); ++j)
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
          }
        } else {
          for (size_t j=0UL; j<result.columns(); ++j) {
            for (size_t i=0UL; i<result.rows(); ++i)
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
          }
        }
        return result ;
      }
    };
     */

    /*

    // Provides only blaze::StaticMatrix<Type,M,N,SO> export
    template<typename Type, size_t M, size_t N, bool SO >
    class Exporter< blaze::StaticMatrix<Type, M, N, SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {
        if ((size_t) mat.nrow() != M && (size_t) mat.ncol() != N) {
          throw std::invalid_argument("Dimension of matrix is not match.");
        }
      }

      blaze::StaticMatrix<Type, M, N, SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        blaze::StaticMatrix<Type, M, N, SO> result;
        if (SO == blaze::rowMajor) {
          for (size_t i=0UL; i<result.rows(); ++i) {
            for( size_t j=0UL; j<result.columns(); ++j) {
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
        } else {
          for (size_t j=0UL; j<result.columns(); ++j) {
            for (size_t i=0UL; i<result.rows(); ++i) {
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
        }
        return result ;
      }
    };

    // Provides only blaze::HybridMatrix<Type,M,N,SO> export
    template<typename Type, size_t M, size_t N, bool SO >
    class Exporter< blaze::HybridMatrix<Type, M, N, SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

    public:
      Exporter(SEXP x) : mat(x) {
        if ((size_t) mat.nrow() != M && (size_t) mat.ncol() != N) {
          throw std::invalid_argument("Dimension of matrix is not match.");
        }
      }

      blaze::HybridMatrix<Type, M, N, SO> get() {
        typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
        blaze::HybridMatrix<Type, M, N, SO> result((size_t) mat.nrow(), (size_t) mat.ncol());
        if (SO == blaze::rowMajor) {
          for (size_t i=0UL; i<result.rows(); ++i) {
            for (size_t j=0UL; j<result.columns(); ++j) {
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
        } else {
          for (size_t j=0UL; j<result.columns(); ++j) {
            for (size_t i=0UL; i<result.rows(); ++i) {
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
        }
        return result ;
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type,blaze::unaligned, blaze::unpadded,SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

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
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

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
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
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
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>( mat(i ,j) );
            }
          }
          return result;
        }
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::unaligned, blaze::padded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type, blaze::unaligned, blaze::padded, SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

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
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>(mat(i ,j));
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
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>(mat(i ,j));
            }
          }
          return result;
        }
      }
    };

    // Provides only blaze::CustomMatrix<Type,blaze::aligned, blaze::padded,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CustomMatrix<Type, blaze::aligned, blaze::padded, SO> > {
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;
      ::Rcpp::Matrix<RTYPE> mat;

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
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>(mat(i ,j));
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
              result(i ,j) = ::Rcpp::internal::caster<value_t, Type>(mat(i ,j));
            }
          }
          return result;
        }
      }
    };

    // Provides only blaze::CompressedMatrix<Type,SO> export
    template<typename Type, bool SO >
    class Exporter< blaze::CompressedMatrix<Type, SO> > {
      ::Rcpp::S4 mat;

    public:
      Exporter(SEXP x) : mat(x) {}
      blaze::CompressedMatrix<Type, SO> get() {
        const int RTYPE = ::Rcpp::traits::r_sexptype_traits<Type>::rtype;

        ::Rcpp::IntegerVector dims = mat.slot("Dim");
        size_t nrow = (size_t) dims[0];
        size_t ncol = (size_t) dims[1];

        // initialize blaze::CompressedMatrix
        blaze::CompressedMatrix<Type,SO> result(nrow, ncol);

        // get the type of sparse matrix
        const char* spMatClass = mat.slot("class");

        // get data of sparse matrix
        ::Rcpp::Vector<RTYPE> x = mat.slot("x");

        if (strcmp(spMatClass + 2, "CMatrix") == 0) {
          ::Rcpp::IntegerVector i = mat.slot("i");
          ::Rcpp::IntegerVector p = mat.slot("p");

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
          ::Rcpp::IntegerVector j = mat.slot("j");
          ::Rcpp::IntegerVector p = mat.slot("p");

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
          ::Rcpp::IntegerVector i = mat.slot("i");
          ::Rcpp::IntegerVector j = mat.slot("j");
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
          ::Rcpp::IntegerVector i = mat.slot("i");
          ::Rcpp::IntegerVector p = mat.slot("p");
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
