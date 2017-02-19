//=================================================================================================
/*!
//  \file blaze/math/lapack/sytrf.h
//  \brief Header file for the LAPACK symmetric matrix decomposition functions (sytrf)
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_LAPACK_SYTRF_H_
#define _BLAZE_MATH_LAPACK_SYTRF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Exception.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/UniqueArray.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
extern "C" {

void ssytrf_( char* uplo, int* n, float*  A, int* lda, int* ipiv, float*  work, int* lwork, int* info );
void dsytrf_( char* uplo, int* n, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info );
void csytrf_( char* uplo, int* n, float*  A, int* lda, int* ipiv, float*  work, int* lwork, int* info );
void zsytrf_( char* uplo, int* n, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK LDLT DECOMPOSITION FUNCTIONS (SYTRF)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LDLT decomposition functions (sytrf) */
//@{
inline void sytrf( char uplo, int n, float* A, int lda, int* ipiv,
                   float* work, int lwork, int* info );

inline void sytrf( char uplo, int n, double* A, int lda, int* ipiv,
                   double* work, int lwork, int* info );

inline void sytrf( char uplo, int n, complex<float>* A, int lda, int* ipiv,
                   complex<float>* work, int lwork, int* info );

inline void sytrf( char uplo, int n, complex<double>* A, int lda, int* ipiv,
                   complex<double>* work, int lwork, int* info );

template< typename MT, bool SO >
inline void sytrf( DenseMatrix<MT,SO>& A, char uplo, int* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the decomposition of the given dense symmetric indefinite single
//        precision column-major matrix.
// \ingroup lapack_decomposition
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix decomposition of a symmetric indefinite single precision
// column-major matrix based on the LAPACK ssytrf() function, which uses the Bunch-Kaufman diagonal
// pivoting method. The decomposition has the form

                      \f[ A = U D U^{T} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but D(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= \a n*NB, where NB
// is the optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a
// workspace query is assumed. The function only calculates the optimal size of the \a work
// array and returns this value as the first entry of the \a work array.
//
// For more information on the ssytrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void sytrf( char uplo, int n, float* A, int lda, int* ipiv,
                   float* work, int lwork, int* info )
{
   ssytrf_( &uplo, &n, A, &lda, ipiv, work, &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the decomposition of the given dense symmetric indefinite double
//        precision column-major matrix.
// \ingroup lapack_decomposition
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix decomposition of a symmetric indefinite double precision
// column-major matrix based on the LAPACK dsytrf() function, which uses the Bunch-Kaufman diagonal
// pivoting method. The decomposition has the form

                      \f[ A = U D U^{T} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but D(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= \a n*NB, where NB
// is the optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a
// workspace query is assumed. The function only calculates the optimal size of the \a work
// array and returns this value as the first entry of the \a work array.
//
// For more information on the dsytrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void sytrf( char uplo, int n, double* A, int lda, int* ipiv,
                   double* work, int lwork, int* info )
{
   dsytrf_( &uplo, &n, A, &lda, ipiv, work, &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the decomposition of the given dense symmetric indefinite single
//        precision complex column-major matrix.
// \ingroup lapack_decomposition
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix decomposition of a symmetric indefinite single precision
// column-major matrix based on the LAPACK csytrf() function, which uses the Bunch-Kaufman diagonal
// pivoting method. The decomposition has the form

                      \f[ A = U D U^{T} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but D(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= \a n*NB, where NB
// is the optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a
// workspace query is assumed. The function only calculates the optimal size of the \a work
// array and returns this value as the first entry of the \a work array.
//
// For more information on the csytrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void sytrf( char uplo, int n, complex<float>* A, int lda, int* ipiv,
                   complex<float>* work, int lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   csytrf_( &uplo, &n, reinterpret_cast<float*>( A ), &lda, ipiv,
            reinterpret_cast<float*>( work ), &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the decomposition of the given dense symmetric indefinite double
//        precision complex column-major matrix.
// \ingroup lapack_decomposition
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix decomposition of a symmetric indefinite double precision
// column-major matrix based on the LAPACK zsytrf() function, which uses the Bunch-Kaufman diagonal
// pivoting method. The decomposition has the form

                      \f[ A = U D U^{T} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but D(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= \a n*NB, where NB
// is the optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a
// workspace query is assumed. The function only calculates the optimal size of the \a work
// array and returns this value as the first entry of the \a work array.
//
// For more information on the zsytrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void sytrf( char uplo, int n, complex<double>* A, int lda, int* ipiv,
                   complex<double>* work, int lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zsytrf_( &uplo, &n, reinterpret_cast<double*>( A ), &lda, ipiv,
            reinterpret_cast<double*>( work ), &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the decomposition of the given dense symmetric indefinite matrix.
// \ingroup lapack_decomposition
//
// \param A The matrix to be decomposed.
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Invalid uplo argument provided.
//
// This function performs the dense matrix decomposition of a symmetric indefinite matrix based
// on the LAPACK sytrf() functions, which use the Bunch-Kaufman diagonal pivoting method. Note
// that the function only works for general, non-adapted matrices with \c float, \c double,
// \c complex<float>, or \c complex<double> element type. The attempt to call the function with
// any adapted matrix or matrices of any other element type results in a compile time error!\n
//
// The decomposition has the form

                      \f[ A = U D U^{T} \texttt{ (if uplo = 'U'), or }
                          A = L D L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U (or \c L) is a product of permutation and unit upper (lower) triangular matrices,
// and \c D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched.
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given \a uplo argument is neither \c 'L' nor \c 'U'.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the sytrf() functions (i.e. ssytrf(), dsytrf(), csytrf(), and zsytrf())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void sytrf( DenseMatrix<MT,SO>& A, char uplo, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   typedef typename MT::ElementType  ET;

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   if( uplo != 'L' && uplo != 'U' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid uplo argument provided" );
   }

   int n   ( numeric_cast<int>( (~A).rows()    ) );
   int lda ( numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   if( n == 0 ) {
      return;
   }

   int lwork( n*lda );
   const UniqueArray<ET> work( new ET[lwork] );

   if( IsRowMajorMatrix<MT>::value ) {
      ( uplo == 'L' )?( uplo = 'U' ):( uplo = 'L' );
   }

   sytrf( uplo, n, (~A).data(), lda, ipiv, work.get(), lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix decomposition" );
}
//*************************************************************************************************

} // namespace blaze

#endif
