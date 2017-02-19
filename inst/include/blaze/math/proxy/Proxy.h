//=================================================================================================
/*!
//  \file blaze/math/proxy/Proxy.h
//  \brief Header file for the Proxy class
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

#ifndef _BLAZE_MATH_PROXY_PROXY_H_
#define _BLAZE_MATH_PROXY_PROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/InversionFlag.h>
#include <blaze/math/proxy/ComplexProxy.h>
#include <blaze/math/proxy/DefaultProxy.h>
#include <blaze/math/proxy/DenseMatrixProxy.h>
#include <blaze/math/proxy/DenseVectorProxy.h>
#include <blaze/math/proxy/SparseMatrixProxy.h>
#include <blaze/math/proxy/SparseVectorProxy.h>
#include <blaze/math/shims/Imaginary.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/shims/IsReal.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Real.h>
#include <blaze/math/traits/AbsExprTrait.h>
#include <blaze/math/traits/AddExprTrait.h>
#include <blaze/math/traits/CTransExprTrait.h>
#include <blaze/math/traits/DivExprTrait.h>
#include <blaze/math/traits/ImagExprTrait.h>
#include <blaze/math/traits/MultExprTrait.h>
#include <blaze/math/traits/RealExprTrait.h>
#include <blaze/math/traits/SubExprTrait.h>
#include <blaze/math/traits/TransExprTrait.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsProxy.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Proxy base class.
// \ingroup math
//
// The Proxy class is a base class for all proxy classes within the \b Blaze library that may
// represent non-numeric data types (vectors, matrices, ...). It augments the interface of the
// deriving proxy class depending on the data type represented by the proxy. In addition, it
// provides an abstraction from the actual type of the proxy, but enables a type-safe conversion
// back to this type via the 'Curiously Recurring Template Pattern' (CRTP).
//
// In order to use the Proxy class it is necessary to publicly derive from it and to provide
// an accessible member function called \a get(), which grants access to the represented element
// via non-const reference. The following example demonstrates these requirements by means of
// the VectorAccessProxy class:

   \code
   template< typename VT >
   class VectorAccessProxy : public Proxy< VectorAccessProxy<VT>, typename VT::ElementType >
   {
      // ...
      typedef typename VT::ElementType  RepresentedType;
      inline RepresentedType& get() const;
      // ...
   };
   \endcode

// The first template parameter specifies the type of the deriving proxy class (CRTP), the second
// template parameter specifies the type of the element represented by the proxy. Within the
// context of the VectorAccessProxy this is the type of the elements of the vector to be accessed.
// Depending on this type the proxy selects the additional interface to provide to the deriving
// class.
*/
template< typename PT          // Type of the proxy
        , typename RT = int >  // Type of the represented element
class Proxy : public If< IsVector<RT>
                       , typename If< IsDenseVector<RT>
                                    , DenseVectorProxy<PT,RT>
                                    , SparseVectorProxy<PT,RT>
                                    >::Type
                       , typename If< IsMatrix<RT>
                                    , typename If< IsDenseMatrix<RT>
                                                 , DenseMatrixProxy<PT,RT>
                                                 , SparseMatrixProxy<PT,RT>
                                                 >::Type
                                    , typename If< IsComplex<RT>
                                                 , ComplexProxy<PT,RT>
                                                 , DefaultProxy<PT,RT>
                                                 >::Type
                                    >::Type
                       >::Type
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Proxy operators */
//@{
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename AddExprTrait<RT1,RT2>::Type
   operator+( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename AddExprTrait<RT,T>::Type >::Type
   operator+( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename AddExprTrait<T,RT>::Type >::Type
   operator+( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename SubExprTrait<RT1,RT2>::Type
   operator-( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename SubExprTrait<RT,T>::Type >::Type
   operator-( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename SubExprTrait<T,RT>::Type >::Type
   operator-( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename MultExprTrait<RT1,RT2>::Type
   operator*( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename MultExprTrait<RT,T>::Type >::Type
   operator*( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename MultExprTrait<T,RT>::Type >::Type
   operator*( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename DivExprTrait<RT1,RT2>::Type
   operator/( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename DivExprTrait<RT,T>::Type >::Type
   operator/( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename DivExprTrait<T,RT>::Type >::Type
   operator/( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator==( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator==( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator==( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator!=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator!=( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator!=( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator<( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator>( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator<=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<=( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<=( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator>=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>=( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>=( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT, typename RT >
inline std::ostream& operator<<( std::ostream& os, const Proxy<PT,RT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the addition.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename AddExprTrait<RT1,RT2>::Type
   operator+( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (~lhs).get() + (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the addition.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename AddExprTrait<RT,T>::Type >::Type
   operator+( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (~lhs).get() + rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the addition.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename AddExprTrait<T,RT>::Type >::Type
   operator+( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs + (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the subtraction.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename SubExprTrait<RT1,RT2>::Type
   operator-( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (~lhs).get() - (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the subtraction.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename SubExprTrait<RT,T>::Type >::Type
   operator-( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (~lhs).get() - rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the subtraction.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename SubExprTrait<T,RT>::Type >::Type
   operator-( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs - (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the multiplication.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename MultExprTrait<RT1,RT2>::Type
   operator*( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (~lhs).get() * (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the multiplication.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename MultExprTrait<RT,T>::Type >::Type
   operator*( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (~lhs).get() * rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the multiplication.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename MultExprTrait<T,RT>::Type >::Type
   operator*( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs * (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the division.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline typename DivExprTrait<RT1,RT2>::Type
   operator/( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (~lhs).get() / (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the division.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, typename DivExprTrait<RT,T>::Type >::Type
   operator/( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (~lhs).get() / rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the division.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, typename DivExprTrait<T,RT>::Type >::Type
   operator/( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs / (~rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator==( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (~lhs).get() == (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator==( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (~lhs).get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator==( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs == (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator!=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (~lhs).get() != (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator!=( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (~lhs).get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator!=( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs != (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator<( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (~lhs).get() < (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (~lhs).get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator>( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (~lhs).get() > (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (~lhs).get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs > (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator<=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (~lhs).get() <= (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<=( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (~lhs).get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator<=( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs <= (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator>=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (~lhs).get() >= (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename PT, typename RT, typename T >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>=( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (~lhs).get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename PT, typename RT >
inline typename DisableIf< IsProxy<T>, bool >::Type
   operator>=( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs >= (~rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the Proxy class template.
// \ingroup math
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant proxy object.
// \return Reference to the output stream.
*/
template< typename PT, typename RT >
inline std::ostream& operator<<( std::ostream& os, const Proxy<PT,RT>& proxy )
{
   return os << (~proxy).get();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Proxy global functions */
//@{
template< typename PT, typename RT >
inline typename TransExprTrait< typename PT::RepresentedType >::Type
   trans( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline typename AbsExprTrait< typename PT::RepresentedType >::Type
   abs( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline typename CTransExprTrait< typename PT::RepresentedType >::Type
   ctrans( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline typename RealExprTrait< typename PT::RepresentedType >::Type
   real( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline typename ImagExprTrait< typename PT::RepresentedType >::Type
   imag( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline void transpose( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline void ctranspose( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline void invert( const Proxy<PT,RT>& proxy );

template< InversionFlag IF, typename PT, typename RT >
inline void invert( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline bool isReal( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline bool isZero( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline bool isOne( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
inline bool isnan( const Proxy<PT,RT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the transpose of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The transpose of the represented element.
//
// This function returns an expression representing the transpose of the element represented by
// the proxy.
*/
template< typename PT, typename RT >
inline typename TransExprTrait< typename PT::RepresentedType >::Typ
    trans( const Proxy<PT,RT>& proxy )
{
   using blaze::trans;

   return trans( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the absolute value of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The absolute value of the represented element.
//
// This function computes the absolute value of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the absolute values of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline typename AbsExprTrait< typename PT::RepresentedType >::Type
   abs( const Proxy<PT,RT>& proxy )
{
   using std::abs;

   return abs( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the conjugate transpose of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The conjugate transpose of the represented element.
//
// This function returns an expression representing the conjugate transpose of the element
// represented by the proxy.
*/
template< typename PT, typename RT >
inline typename CTransExprTrait< typename PT::RepresentedType >::Type
   ctrans( const Proxy<PT,RT>& proxy )
{
   using blaze::ctrans;

   return ctrans( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the real part of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The real part of the represented element.
//
// This function returns the real part of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the real part of each each element of the vector/matrix.
*/
template< typename PT, typename RT >
inline typename RealExprTrait< typename PT::RepresentedType >::Type
   real( const Proxy<PT,RT>& proxy )
{
   using blaze::real;

   return real( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the imaginary part of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The imaginary part of the represented element.
//
// This function returns the imaginary part of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the real part of each each element of the vector/matrix.
*/
template< typename PT, typename RT >
inline typename ImagExprTrait< typename PT::RepresentedType >::Type
   imag( const Proxy<PT,RT>& proxy )
{
   using blaze::imag;

   return imag( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place transpose of the represented matrix element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::logic_error Matrix cannot be transposed.
//
// This function transposes the represented matrix in-place. The transpose operation fails if ...
//
//  - ... the represented matrix has a fixed size and is non-square;
//  - ... the represented matrix is a triangular matrix.
//
// In all failure cases a \a std::logic_error exception is thrown. Additionally, in case the
// represented matrix cannot be modified, a \a std::invalid_argument exception is thrown.
*/
template< typename PT, typename RT >
inline void transpose( const Proxy<PT,RT>& proxy )
{
   if( (~proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   transpose( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place conjugate transpose of the represented matrix element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::logic_error Matrix cannot be transposed.
//
// This function transposes the represented matrix in-place. The transpose operation fails if ...
//
//  - ... the represented matrix has a fixed size and is non-square;
//  - ... the represented matrix is a triangular matrix.
//
// In all failure cases a \a std::logic_error exception is thrown. Additionally, in case the
// represented matrix cannot be modified, a \a std::invalid_argument exception is thrown.
*/
template< typename PT, typename RT >
inline void ctranspose( const Proxy<PT,RT>& proxy )
{
   if( (~proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   ctranspose( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the represented scalar or dense matrix element. The inversion fails if
// the represented element is a dense matrix, which ...
//
//  - ... is not a square matrix;
//  - ... is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \a std::invalid_argument exception is thrown. Additionally, in case the
// represented scalar or matrix cannot be modified, a \a std::invalid_argument exception is thrown.
//
// \note In case the represented element is a dense matrix, this function does not provide any
// exception safety guarantee, i.e. in case an exception is thrown the matrix may already have
// been modified.
//
// \note In case the represented element is a dense matrix, this function can only be used if the
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error will
// be created.
*/
template< typename PT, typename RT >
inline void invert( const Proxy<PT,RT>& proxy )
{
   if( (~proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   invert( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the represented dense matrix element by means of the specified matrix
// inversion algorithm \c IF:

   \code
   invert<byLU>( A );    // Inversion of a general matrix
   invert<byLDLT>( A );  // Inversion of a symmetric indefinite matrix
   invert<byLDLH>( A );  // Inversion of a Hermitian indefinite matrix
   invert<byLLH>( A );   // Inversion of a Hermitian positive definite matrix
   \endcode

// The inversion fails if the represented dense matrix element ...
//
//  - ... is not a square matrix;
//  - ... is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \a std::invalid_argument exception is thrown. Additionally, in case the
// represented scalar or matrix cannot be modified, a \a std::invalid_argument exception is thrown.
//
// \note In case the represented element is a dense matrix, this function does not provide any
// exception safety guarantee, i.e. in case an exception is thrown the matrix may already have
// been modified.
//
// \note In case the represented element is a dense matrix, this function can only be used if the
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error will
// be created.
*/
template< InversionFlag IF, typename PT, typename RT >
inline void invert( const Proxy<PT,RT>& proxy )
{
   if( (~proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   invert<IF>( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the element represents a real number.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the element represents a real number, \a false otherwise.
//
// This function checks whether the element represented by the proxy represents the a real
// number. In case the element is of built-in type, the function returns \a true. In case
// the element is of complex type, the function returns \a true if the imaginary part is
// equal to 0. Otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isReal( const Proxy<PT,RT>& proxy )
{
   using blaze::isReal;

   return isReal( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is 0.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is 0, \a false otherwise.
//
// This function checks whether the element represented by the proxy represents the numeric
// value 0. In case it is 0, the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isZero( const Proxy<PT,RT>& proxy )
{
   using blaze::isZero;

   return isZero( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is 1.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is 1, \a false otherwise.
//
// This function checks whether the element represented by the proxy represents the numeric
// value 1. In case it is 1, the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isOne( const Proxy<PT,RT>& proxy )
{
   using blaze::isOne;

   return isOne( (~proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is not a number.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is in not a number, \a false otherwise.
//
// This function checks whether the element represented by the proxy is not a number (NaN).
// In case it is not a number, the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isnan( const Proxy<PT,RT>& proxy )
{
   using blaze::isnan;

   return isnan( (~proxy).get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
