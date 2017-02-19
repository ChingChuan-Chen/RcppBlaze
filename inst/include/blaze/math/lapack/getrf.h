//=================================================================================================
/*!
//  \file blaze/math/lapack/getrf.h
//  \brief Header file for the LAPACK LU decomposition functions (getrf)
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

#ifndef _BLAZE_MATH_LAPACK_GETRF_H_
#define _BLAZE_MATH_LAPACK_GETRF_H_


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
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
extern "C" {

void sgetrf_( int* m, int* n, float*  A, int* lda, int* ipiv, int* info );
void dgetrf_( int* m, int* n, double* A, int* lda, int* ipiv, int* info );
void cgetrf_( int* m, int* n, float*  A, int* lda, int* ipiv, int* info );
void zgetrf_( int* m, int* n, double* A, int* lda, int* ipiv, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK LU DECOMPOSITION FUNCTIONS (GETRF)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LU decomposition functions (getrf) */
//@{
inline void getrf( int m, int n, float* A, int lda, int* ipiv, int* info );

inline void getrf( int m, int n, double* A, int lda, int* ipiv, int* info );

inline void getrf( int m, int n, complex<float>* A, int lda, int* ipiv, int* info );

inline void getrf( int m, int n, complex<double>* A, int lda, int* ipiv, int* info );

template< typename MT, bool SO >
inline void getrf( DenseMatrix<MT,SO>& A, int* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general single precision
//        column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general m-by-n single precision
// column-major matrix based on the LAPACK sgetrf() function, which uses partial pivoting with row
// interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the sgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getrf( int m, int n, float* A, int lda, int* ipiv, int* info )
{
   sgetrf_( &m, &n, A, &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general double precision
//        column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general m-by-n double precision
// column-major matrix based on the LAPACK dgetrf() function, which uses partial pivoting with row
// interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the dgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getrf( int m, int n, double* A, int lda, int* ipiv, int* info )
{
   dgetrf_( &m, &n, A, &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general single precision
//        complex column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general m-by-n single precision
// complex column-major matrix based on the LAPACK cgetrf() function, which uses partial pivoting
// with row interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the cgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getrf( int m, int n, complex<float>* A, int lda, int* ipiv, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cgetrf_( &m, &n, reinterpret_cast<float*>( A ), &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general double precision
//        complex column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general m-by-n double precision
// complex column-major matrix based on the LAPACK zgetrf() function, which uses partial pivoting
// with row interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the zgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getrf( int m, int n, complex<double>* A, int lda, int* ipiv, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zgetrf_( &m, &n, reinterpret_cast<double*>( A ), &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief LAPACK kernel for the LU decomposition of the given dense general matrix.
// \ingroup lapack_decomposition
//
// \param A The matrix to be decomposed.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \return void
//
// This function performs the dense matrix LU decomposition of a general m-by-n matrix based
// on the LAPACK \c getrf() functions, which use partial pivoting with row/column interchanges.
// Note that the function only works for general, non-adapted matrices with \c float, \c double,
// \c complex<float>, or \c complex<double> element type. The attempt to call the function with
// adaptors or matrices of any other element type results in a compile time error!\n
//
// In case of a column-major matrix, the resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c L is a lower unitriangular matrix (lower trapezoidal if \a m > \a n), \c U is an upper
// triangular matrix (upper trapezoidal if \a m < \a n), and \c P is an m-by-m permutation matrix,
// which represents the pivoting indices for the applied row interchanges.
//
// In case of a row-major matrix, the resulting decomposition has the form

                          \f[ A = L \cdot U \cdot P, \f]

// where \c P is an n-by-n permutation matrix, which represents the pivoting indices for the applied
// column interchanges, \c L is a lower triangular matrix (lower trapezoidal if \a m > \a n), and
// \c U is an upper unitriangular matrix (upper trapezoidal if \a m < \a n).
//
// The resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part
// of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L or \c U are
// not stored.
//
// For more information on the getrf() functions (i.e. sgetrf(), dgetrf(), cgetrf(), and zgetrf())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
//
// \note The LU decomposition will never fail, even for singular matrices. However, in case of a
// singular matrix the resulting decomposition cannot be used for a matrix inversion or solving
// a linear system of equations.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void getrf( DenseMatrix<MT,SO>& A, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   int m   ( numeric_cast<int>( SO ? (~A).rows() : (~A).columns() ) );
   int n   ( numeric_cast<int>( SO ? (~A).columns() : (~A).rows() ) );
   int lda ( numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   if( m == 0 || n == 0 ) {
      return;
   }

   getrf( m, n, (~A).data(), lda, ipiv, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for LU decomposition" );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
