//=================================================================================================
/*!
//  \file blaze/math/lapack/orgrq.h
//  \brief Header file for the LAPACK functions to reconstruct Q from a RQ decomposition (orgrq)
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

#ifndef _BLAZE_MATH_LAPACK_ORGRQ_H_
#define _BLAZE_MATH_LAPACK_ORGRQ_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Builtin.h>
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

void sorgrq_( int* m, int* n, int* k, float*  A, int* lda, float*  tau, float*  work, int* lwork, int* info );
void dorgrq_( int* m, int* n, int* k, double* A, int* lda, double* tau, double* work, int* lwork, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK FUNCTIONS TO RECONSTRUCT Q FROM A RQ DECOMPOSITION (ORGRQ)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK functions to reconstruct Q from a RQ decomposition (orgrq) */
//@{
inline void orgrq( int m, int n, int k, float* A, int lda, const float* tau,
                   float* work, int lwork, int* info );

inline void orgrq( int m, int n, int k, double* A, int lda, const double* tau,
                   double* work, int lwork, int* info );

template< typename MT, bool SO >
inline void orgrq( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the reconstruction of the orthogonal matrix Q from a RQ decomposition.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..n)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param k The number of elementary reflectors, whose product defines the matrix \f$[0..m)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function generates all or part of the orthogonal matrix Q from a RQ decomposition based on
// the LAPACK sorgrq() function for single precision column-major matrices that have already been
// factorized by the sgerqf() function. The \a info argument provides feedback on the success of
// the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the sorgrq() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void orgrq( int m, int n, int k, float* A, int lda, const float* tau, float* work, int lwork, int* info )
{
   sorgrq_( &m, &n, &k, A, &lda, const_cast<float*>( tau ), work, &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the reconstruction of the orthogonal matrix Q from a RQ decomposition.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..n)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param k The number of elementary reflectors, whose product defines the matrix \f$[0..m)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function generates all or part of the orthogonal matrix Q from a RQ decomposition based on
// the LAPACK dorgrq() function for double precision column-major matrices that have already been
// factorized by the dgerqf() function. The \a info argument provides feedback on the success of
// the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the dorgrq() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void orgrq( int m, int n, int k, double* A, int lda, const double* tau, double* work, int lwork, int* info )
{
   dorgrq_( &m, &n, &k, A, &lda, const_cast<double*>( tau ), work, &lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the reconstruction of the orthogonal matrix Q from a RQ decomposition.
// \ingroup lapack_decomposition
//
// \param A The decomposed matrix.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \return void
//
// This function reconstructs the orthogonal matrix \a Q of a RQ decomposition based on the LAPACK
// orgrq() functions from matrices that have already been RQ factorized by the gerqf() functions.
// Note that this function can only be used for general, non-adapted matrices with \c float or
// \c double element type. The attempt to call the function with any adapted matrix or matrices
// of any other element type results in a compile time error!\n
//
// The row-major \a m-by-min(\a m,\a n) or column-major min(\a m,\a n)-by-\a n \a Q matrix is
// stored in the within the given matrix \a A:

   \code
   using blaze::DynamicMatrix;
   using blaze::columnMajor;

   DynamicMatrix<double,columnMajor> A;
   DynamicVector<double> tau;
   // ... Resizing and initialization

   gerqf( A, tau.data() );  // Performing the RQ decomposition
   orgrq( A, tau.data() );  // Reconstructing the Q matrix

   const int m( A.rows() );
   const int n( A.columns() );

   const size_t row( m > n ? m - n : 0UL )
   DynamicMatrix<double,columnMajor> Q( submatrix( A, row, 0UL, min(m,n), n ) );
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::rowMajor;

   DynamicMatrix<double,rowMajor> A;
   DynamicVector<double> tau;
   // ... Resizing and initialization

   gerqf( A, tau.data() );  // Performing the RQ decomposition
   orgrq( A, tau.data() );  // Reconstructing the Q matrix

   const int m( A.rows() );
   const int n( A.columns() );

   const size_t column( m < n ? n - m : 0UL )
   DynamicMatrix<double,rowMajor> Q( submatrix( A, 0UL, column, m, min(m,n) ) );
   \endcode

// For more information on the orgrq() functions (i.e. sorgrq() and dorgrq()) see the LAPACK
// online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT, bool SO >
inline void orgrq( DenseMatrix<MT,SO>& A, const typename MT::ElementType* tau )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( typename MT::ElementType );

   typedef typename MT::ElementType  ET;

   int n   ( numeric_cast<int>( SO ? (~A).columns() : (~A).rows() ) );
   int m   ( numeric_cast<int>( SO ? (~A).rows() : (~A).columns() ) );
   int k   ( min( m, n ) );
   int lda ( numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   if( k == 0 ) {
      return;
   }

   int lwork( k*lda );
   const UniqueArray<ET> work( new ET[lwork] );
   const size_t offset( ( m > n )?( m - n ):( 0UL ) );

   orgrq( k, n, k, (~A).data()+offset, lda, tau, work.get(), lwork, &info );

   BLAZE_INTERNAL_ASSERT( info == 0, "Invalid argument for Q reconstruction" );
}
//*************************************************************************************************

} // namespace blaze

#endif
