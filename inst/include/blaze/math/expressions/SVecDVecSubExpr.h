//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecDVecSubExpr.h
//  \brief Header file for the sparse vector/dense vector subtraction expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECDVECSUBEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECDVECSUBEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/constraints/VecVecSubExpr.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecVecSubExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/Max.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SVECDVECSUBEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse vector-dense vector subtractions.
// \ingroup dense_vector_expression
//
// The SVecDVecSubExpr class represents the compile time expression for subtractions between
// a sparse vector and a dense vector.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
class SVecDVecSubExpr : public DenseVector< SVecDVecSubExpr<VT1,VT2,TF>, TF >
                      , private VecVecSubExpr
                      , private Computation
{
 private:
   //**Type definitions****************************************************************************
   typedef typename VT1::ResultType     RT1;  //!< Result type of the left-hand side sparse vector expression.
   typedef typename VT2::ResultType     RT2;  //!< Result type of the right-hand side dense vector expression.
   typedef typename VT1::ReturnType     RN1;  //!< Return type of the left-hand side sparse vector expression.
   typedef typename VT2::ReturnType     RN2;  //!< Return type of the right-hand side dense vector expression.
   typedef typename VT1::CompositeType  CT1;  //!< Composite type of the left-hand side sparse vector expression.
   typedef typename VT2::CompositeType  CT2;  //!< Composite type of the right-hand side dense vector expression.
   typedef typename VT1::TransposeType  TT1;  //!< Transpose type of the left-hand side sparse vector expression.
   typedef typename VT2::TransposeType  TT2;  //!< Transpose type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   enum { returnExpr = !IsTemporary<RN1>::value && !IsTemporary<RN2>::value };

   //! Expression return type for the subscript operator.
   typedef typename SubExprTrait<RN1,RN2>::Type  ExprReturnType;
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the explicit application of the SFINAE principle.
   /*! The UseSMPAssign struct is a helper struct for the selection of the parallel evaluation
       strategy. In case at least one of the two vector operands is not SMP assignable, \a value
       is set to 1 and the expression specific evaluation strategy is selected. Otherwise \a value
       is set to 0 and the default strategy is chosen. */
   template< typename VT >
   struct UseSMPAssign {
      enum { value = ( !VT1::smpAssignable || !VT2::smpAssignable ) };
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   typedef SVecDVecSubExpr<VT1,VT2,TF>         This;           //!< Type of this SVecDVecSubExpr instance.
   typedef typename SubTrait<RT1,RT2>::Type    ResultType;     //!< Result type for expression template evaluations.
   typedef typename ResultType::TransposeType  TransposeType;  //!< Transpose type for expression template evaluations.
   typedef typename ResultType::ElementType    ElementType;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   typedef const typename SelectType< returnExpr, ExprReturnType, ElementType >::Type  ReturnType;

   //! Data type for composite expression templates.
   typedef const ResultType  CompositeType;

   //! Composite type of the left-hand side sparse vector expression.
   typedef typename SelectType< IsExpression<VT1>::value, const VT1, const VT1& >::Type  LeftOperand;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename SelectType< IsExpression<VT2>::value, const VT2, const VT2& >::Type  RightOperand;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   enum { vectorizable = 0 };

   //! Compilation switch for the expression template assignment strategy.
   enum { smpAssignable = 0 };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SVecDVecSubExpr class.
   //
   // \param lhs The left-hand side sparse vector operand of the subtraction expression.
   // \param rhs The right-hand side dense vector operand of the subtraction expression.
   */
   explicit inline SVecDVecSubExpr( const VT1& lhs, const VT2& rhs )
      : lhs_( lhs )  // Left-hand side sparse vector of the subtraction expression
      , rhs_( rhs )  // Right-hand side dense vector of the subtraction expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.size() == rhs.size(), "Invalid vector sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < lhs_.size(), "Invalid vector access index" );
      return lhs_[index] - rhs_[index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= lhs_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const {
      return lhs_.size();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse vector operand.
   //
   // \return The left-hand side sparse vector operand.
   */
   inline LeftOperand leftOperand() const {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const {
      return rhs_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const {
      return ( lhs_.canAlias( alias ) ) ||
             ( IsExpression<VT2>::value && rhs_.canAlias( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side sparse vector of the subtraction expression.
   RightOperand rhs_;  //!< Right-hand side dense vector of the subtraction expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-dense
   // vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign   ( ~lhs, -rhs.rhs_ );
      addAssign( ~lhs,  rhs.lhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-dense vector addition to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-dense
   // vector addition expression to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse vector-
   // dense vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( ~lhs, rhs.rhs_ );
      addAssign( ~lhs, rhs.lhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse vector-
   // dense vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( ~lhs, rhs.rhs_ );
      subAssign( ~lhs, rhs.lhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // vector-dense vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse vector-dense
   // vector subtraction expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT> >::Type
      smpAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign   ( ~lhs, -rhs.rhs_ );
      smpAddAssign( ~lhs,  rhs.lhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse vector-dense vector addition to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse vector-dense
   // vector addition expression to a sparse vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline typename EnableIf< UseSMPAssign<VT> >::Type
      smpAssign( SparseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // vector-dense vector subtraction expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT> >::Type
      smpAddAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpSubAssign( ~lhs, rhs.rhs_ );
      smpAddAssign( ~lhs, rhs.lhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse vector-dense vector subtraction to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // vector-dense vector subtraction expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT> >::Type
      smpSubAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAddAssign( ~lhs, rhs.rhs_ );
      smpSubAssign( ~lhs, rhs.lhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a sparse vector-dense vector subtraction to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // sparse vector-dense vector subtraction expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline typename EnableIf< UseSMPAssign<VT> >::Type
      smpMultAssign( DenseVector<VT,TF>& lhs, const SVecDVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( typename ResultType::CompositeType );

      BLAZE_INTERNAL_ASSERT( (~lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT1, TF );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_VECVECSUBEXPR( VT1, VT2 );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of a sparse vector and a dense vector
//        (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector for the vector subtraction.
// \param rhs The right-hand side dense vector to be subtracted from the sparse vector.
// \return The difference of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the subtraction of a sparse vector and a dense vector:

   \code
   blaze::CompressedVector<double> a;
   blaze::DynamicVector<double> b, c;
   // ... Resizing and initialization
   c = a - b;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved vector element types \a T1::ElementType and \a T2::ElementType.
// Both vector types \a T1 and \a T2 as well as the two element types \a T1::ElementType and
// \a T2::ElementType have to be supported by the SubTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1  // Type of the left-hand side sparse vector
        , typename T2  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag
inline const SVecDVecSubExpr<T1,T2,TF>
   operator-( const SparseVector<T1,TF>& lhs, const DenseVector<T2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).size() != (~rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   return SVecDVecSubExpr<T1,T2,TF>( ~lhs, ~rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition operator for the addition of a sparse vector-dense vector subtraction
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}-\vec{c})+\vec{d} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector-dense vector subtraction.
// \param rhs The right-hand side dense vector.
// \return The sum of the two vectors.
//
// This operator implements a performance optimized treatment of the addition of a sparse
// vector-dense vector subtraction expression to a dense vector.
*/
template< typename T1    // Type of the sparse vector of the left-hand side expression
        , typename T2    // Type of the dense vector of the left-hand side expression
        , bool TF        // Transpose flag of the left-hand side expression
        , typename T3 >  // Type of right-hand side dense vector
inline const typename AddExprTrait< SVecDVecSubExpr<T1,T2,TF>, T3 >::Type
   operator+( const SVecDVecSubExpr<T1,T2,TF>& lhs, const DenseVector<T3,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (~rhs) - lhs.rightOperand() ) + lhs.leftOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction operator for the subtraction of a sparse vector-dense vector subtraction
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}-\vec{c})-\vec{d} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector-dense vector subtraction.
// \param rhs The right-hand side dense vector.
// \return The difference of the two vectors.
//
// This operator implements a performance optimized treatment of the subtraction of a
// sparse vector-dense vector subtraction expression and a dense vector.
*/
template< typename T1    // Type of the sparse vector of the left-hand side expression
        , typename T2    // Type of the dense vector of the left-hand side expression
        , bool TF        // Transpose flag of the left-hand side expression
        , typename T3 >  // Type of right-hand side dense vector
inline const typename SubExprTrait< SVecDVecSubExpr<T1,T2,TF>, T3 >::Type
   operator-( const SVecDVecSubExpr<T1,T2,TF>& lhs, const DenseVector<T3,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return lhs.leftOperand() - ( lhs.rightOperand() + (~rhs) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool TF >
struct Size< SVecDVecSubExpr<VT1,VT2,TF> >
   : public Max< Size<VT1>, Size<VT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct DVecDVecAddExprTrait< SVecDVecSubExpr<VT1,VT2,false>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsSparseVector<VT1>::value && IsColumnVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsColumnVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && IsColumnVector<VT3>::value
                              , typename DVecSVecAddExprTrait< typename DVecDVecSubExprTrait<VT3,VT2>::Type, VT1 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct TDVecTDVecAddExprTrait< SVecDVecSubExpr<VT1,VT2,true>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsSparseVector<VT1>::value && IsRowVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsRowVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && IsRowVector<VT3>::value
                              , typename TDVecTSVecAddExprTrait< typename TDVecTDVecSubExprTrait<VT3,VT2>::Type, VT1 >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct DVecDVecSubExprTrait< SVecDVecSubExpr<VT1,VT2,false>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsSparseVector<VT1>::value && IsColumnVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsColumnVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && IsColumnVector<VT3>::value
                              , typename SVecDVecSubExprTrait< VT1, typename DVecDVecAddExprTrait<VT2,VT3>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename VT3 >
struct TDVecTDVecSubExprTrait< SVecDVecSubExpr<VT1,VT2,true>, VT3 >
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename SelectType< IsSparseVector<VT1>::value && IsRowVector<VT1>::value &&
                                IsDenseVector<VT2>::value  && IsRowVector<VT2>::value &&
                                IsDenseVector<VT3>::value  && IsRowVector<VT3>::value
                              , typename TSVecTDVecSubExprTrait< VT1, typename TDVecTDVecAddExprTrait<VT2,VT3>::Type >::Type
                              , INVALID_TYPE >::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool TF, bool AF >
struct SubvectorExprTrait< SVecDVecSubExpr<VT1,VT2,TF>, AF >
{
 public:
   //**********************************************************************************************
   typedef typename SubExprTrait< typename SubvectorExprTrait<const VT1,AF>::Type
                                , typename SubvectorExprTrait<const VT2,AF>::Type >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
