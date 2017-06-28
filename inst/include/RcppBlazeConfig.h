// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
// Copyright (C)  2017         Chingchuan Chen
//
// This file is based on RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze.
//
// RcppBlazeConfig.h: Rcpp/Blaze glue
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
#pragma GCC diagnostic ignored "-Wpedantic"

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
