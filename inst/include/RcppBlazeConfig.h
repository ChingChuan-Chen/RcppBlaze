// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// RcppBlazeConfig.h: Rcpp/Blaze glue
//
// Copyright (C)  2017 Chingchuan Chen
//
// This file is part of RcppBlaze.
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

#ifndef RcppBlaze__RcppBlazeConfig__h
#define RcppBlaze__RcppBlazeConfig__h

// remove some warnings
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#define BOOST_ERROR_CODE_HEADER_ONLY
#define BOOST_SYSTEM_NO_DEPRECATED

// avoid conflict to R constants
#define _BLAZE_MATH_CONSTANTS_H_

#if defined(BLAZE_USE_BOOST_THREADS)
#error "Boost threads could not be used!"
#endif

#if defined(_OPENMP) && defined(BLAZE_USE_CPP_THREADS)
#error "Openmp is used, not to use C++11 threads!"
#undef BLAZE_USE_CPP_THREADS
#endif

#endif
