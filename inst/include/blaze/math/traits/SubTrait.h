//=================================================================================================
/*!
//  \file blaze/math/traits/SubTrait.h
//  \brief Header file for the subtraction trait
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

#ifndef _BLAZE_MATH_TRAITS_SUBTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/typeof/typeof.hpp>
#include <blaze/util/Complex.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/CommonType.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the SubTrait class.
// \ingroup math_traits
//
// \section subtrait_general General
//
// The SubTrait class template offers the possibility to select the resulting data type of a
// generic subtraction operation between the two given types \a T1 and \a T2. SubTrait defines
// the nested type \a Type, which represents the resulting data type of the subtraction. In case
// the two types \a T1 and \a T2 cannot be subtracted, a compilation error is created. Note that
// \c const and \c volatile qualifiers and reference modifiers are generally ignored.
//
// Per default, SubTrait supports all built-in data types. Additionally, the Blaze library
// provides appropriate specializations for the following user-defined arithmetic types:
//
// <ul>
//    <li>std::complex</li>
//    <li>blaze::StaticVector</li>
//    <li>blaze::HybridVector</li>
//    <li>blaze::DynamicVector</li>
//    <li>blaze::CustomVector</li>
//    <li>blaze::CompressedVector</li>
//    <li>blaze::StaticMatrix</li>
//    <li>blaze::HybridMatrix</li>
//    <li>blaze::DynamicMatrix</li>
//    <li>blaze::CustomMatrix</li>
//    <li>blaze::CompressedMatrix</li>
//    <li>blaze::SymmetricMatrix</li>
//    <li>blaze::HermitianMatrix</li>
//    <li>blaze::LowerMatrix</li>
//    <li>blaze::UniLowerMatrix</li>
//    <li>blaze::StrictlyLowerMatrix</li>
//    <li>blaze::UpperMatrix</li>
//    <li>blaze::UniUpperMatrix</li>
//    <li>blaze::StrictlyUpperMatrix</li>
//    <li>blaze::DiagonalMatrix</li>
// </ul>
//
//
// \n \section subtrait_specializations Creating custom specializations
//
// AddTrait is guaranteed to work for all data types that provide a subtraction operator (i.e.
// \c operator-). In order to add support for user-defined data types that either don't provide
// a subtraction operator or whose subtraction operator returns a proxy object instead of a
// concrete type (as it is for instance common in expression template libraries) it is possible
// to specialize the SubTrait template. The following example shows the according specialization
// for the subtraction between two dynamic column vectors:

   \code
   template< typename T1, typename T2 >
   struct SubTrait< DynamicVector<T1,columnVector>, DynamicVector<T2,columnVector> >
   {
      typedef DynamicVector< typename SubTrait<T1,T2>::Type, columnVector >  Type;
   };
   \endcode

// \n \section subtrait_examples Examples
//
// The following example demonstrates the use of the SubTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >  // The two generic types
   typename SubTrait<T1,T2>::Type        // The resulting generic return type
   sub( const T1& t1, const T2& t2 )     //
   {                                     // The function 'sub' returns the
      return t1 - t2;                    // difference of the two given values
   }                                     //
   \endcode
*/
template< typename T1        // Type of the left-hand side operand
        , typename T2        // Type of the right-hand side operand
        , typename = void >  // Restricting condition
struct SubTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename RemoveReference< typename RemoveCV<T1>::Type >::Type  Type1;
   typedef typename RemoveReference< typename RemoveCV<T2>::Type >::Type  Type2;
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct SubType { typedef BOOST_TYPEOF_TPL( Type1() - Type2() )  Type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename If< Or< IsConst<T1>, IsVolatile<T1>, IsReference<T1>
                          , IsConst<T2>, IsVolatile<T2>, IsReference<T2> >
                      , SubTrait<Type1,Type2>, SubType >::Type::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SubTrait class template for a complex and a built-in type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct SubTrait< complex<T1>, T2, typename EnableIf< IsBuiltin<T2> >::Type >
{
 public:
   //**********************************************************************************************
   typedef typename CommonType< complex<T1> , T2 >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SubTrait class template for a built-in and a complex type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct SubTrait< T1, complex<T2>, typename EnableIf< IsBuiltin<T1> >::Type >
{
 public:
   //**********************************************************************************************
   typedef typename CommonType< T1, complex<T2> >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SubTrait class template for two complex types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct SubTrait< complex<T1>, complex<T2> >
{
 public:
   //**********************************************************************************************
   typedef typename CommonType< complex<T1>, complex<T2> >::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
